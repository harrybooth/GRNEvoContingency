using DifferentialEquations
using Distributions
using Distributed
using DiffEqBase
using StatsBase
using Random

# Solvers 

struct DESystemSolver{A <: DEAlgorithm}
    alg :: A
    kwargs :: NamedTuple
end

function DefaultGRNSolver()
    DESystemSolver(Tsit5(),(isoutofdomain=(u,p,t) -> any(x -> x < 0, u), reltol = 1e-6,abstol = 1e-8,callback = TerminateSteadyState(1e-8,1e-6),maxiters = 1e3, verbose = false, save_everystep = false))
end

function TimeStampedGRNSolver(save_at)
    DESystemSolver(Tsit5(),(isoutofdomain=(u,p,t) -> any(x -> x < 0, u), reltol = 1e-6,abstol = 1e-8,maxiters = 1e3, verbose = false, saveat = save_at))
end

function TimeStampedGRNSolver(save_at,save_id)
    DESystemSolver(Tsit5(),(isoutofdomain=(u,p,t) -> any(x -> x < 0, u), reltol = 1e-6,abstol = 1e-8,maxiters = 1e3, verbose = false, saveat = save_at,save_idxs = save_id))
end

function DenseGRNSolver(save_id)
    DESystemSolver(Tsit5(),(isoutofdomain=(u,p,t) -> any(x -> x < 0, u), reltol = 1e-6,abstol = 1e-8,maxiters = 1e3, verbose = false, dense = true,save_idxs = save_id))
end

# Model parameters

struct GRNParameters
    degradation :: Vector{Float64}
    g0 :: Matrix{Float64}
end

function DefaultGRNParameters()
    GRNParameters(deg_rate_g .* ones(Ng),init_conc_g .* ones(Ng,Nc))
end

function NUGRNParameters(id_large::Vector{Int})
    g_init = ones(Ng,Nc)

    id_small = [i for i in 1:3 if i ∉ id_large]

    for i in id_small
        g_init[i,:] .= g_init[i,:] .* init_conc_g
    end

    for i in id_large
        g_init[i,:] .= g_init[i,:] .* init_conc_g_large 
    end

    GRNParameters(deg_rate_g .* ones(Ng),g_init)
end

# Individual and populations

struct Individual
    genotype :: DEProblem
    phenotype :: DESolution
end

function Individual(genotype::DEProblem,development::DESystemSolver)
    phenotype  = solve(genotype,development.alg;development.kwargs...)
    Individual(genotype,phenotype)
end

function Individual(start_network::Matrix{Float64},grn_parameters::GRNParameters,development::DESystemSolver)

    p = (start_network,grn_parameters.degradation)

    genotype = ODEProblem(gene_regulation_1d!,grn_parameters.g0,(0,Inf),p)
    phenotype  = solve(genotype,development.alg;development.kwargs...)
    
    Individual(genotype,phenotype)
end

function Individual(start_network::Matrix{Float64},t2s::Float64,grn_parameters::GRNParameters,development::DESystemSolver)

    p = (start_network,grn_parameters.degradation)

    genotype = ODEProblem(gene_regulation_1d!,grn_parameters.g0,(0,t2s+eps()),p)
    phenotype  = solve(genotype,development.alg;development.kwargs...)
    
    Individual(genotype,phenotype)
end

mutable struct Population{T}
    dominant_individual::Individual
    fitness :: T
    has_fixed :: Bool
end

# Mutation

struct MutationOperatorDual
    mult_noise_distribution :: Distribution
    additive_noise_distribution :: Distribution
    n_sample_func :: Any
    pm_prob :: Float64
    start_affinity :: Float64
    max_w ::Float64
    mutation_weights :: Vector{CartesianIndex{2}}
    sign_flip_probability :: Float64
end

function create_mutant(ind::Individual,mutate_function,development)
    new_w, m_choices,m_type,m_sizes,valid = mutate_function(ind.genotype.p[1])
    Individual(remake(ind.genotype, p = (new_w,ind.genotype.p[2:end]...)),development),m_choices,m_type,m_sizes,valid
end

function create_mutant(ind::Individual,new_w::Matrix{Float64},development)
    Individual(remake(ind.genotype, p = (new_w,ind.genotype.p[2:end]...)),development),nothing,nothing,nothing,nothing
end

function noise_mtype_mult_add(w::Matrix{Float64},mut_op::MutationOperatorDual)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    mut_op

    choices = sort(StatsBase.sample(mut_op.mutation_weights,n_mut,replace = false))
    mtype = []
    sizes = []

    for index in choices

        if rand() < mut_op.pm_prob

            push!(mtype,:multiplicative)

            if new_w[index] == 0
                n = rand(mut_op.mult_noise_distribution)
                if rand() < 0.5
                    new_w[index] = mut_op.start_affinity*n
                    push!(sizes,n)
                else
                    new_w[index] = -1*mut_op.start_affinity*n
                    push!(sizes,-n)
                end
            else
                n = rand(mut_op.mult_noise_distribution)
                new_w[index] = new_w[index]*n
                push!(sizes,n)
            end

        else
            push!(mtype,:additive)

            n = rand(mut_op.additive_noise_distribution)

            new_w[index] = new_w[index] + n
            push!(sizes,n)
        end

        if abs(new_w[index]) > mut_op.max_w
            new_w[index] = sign(new_w[index])*mut_op.max_w
        end

    end

    return new_w, choices, mtype, sizes, true
end

# Selection 

function fixation_probability_kim(Δf,β,N)
    Δf != 0 ? (1 - exp(-2*β*Δf)) / (1 - exp(-2*β*N*Δf)) : 1/N
end

function fixation_probability_kim(Δf1,Δf2,β,N)
    Δf1 != 0 ? (1 - exp(-2*β*Δf1)) / (1 - exp(-2*β*N*Δf1)) : Δf2 != 0 ? (1 - exp(-2*β*Δf2)) / (1 - exp(-2*β*N*Δf2)) : 1/N
end

function strong_selection!(population::Population{Tuple{Float64,Float64}},mutant::Individual,β::Tuple{Float64,Int64},fitness_function)

    mutant_fitness = fitness_function(mutant.phenotype)

    population.has_fixed = false

    if rand() < fixation_probability_kim(mutant_fitness[1] - population.fitness[1],mutant_fitness[2] - population.fitness[2],β[1],β[2])
        population.dominant_individual = mutant
        population.fitness = mutant_fitness
        population.has_fixed = true
    end
end

# Fitness fitness_evaluation

function fitness_evaluation(sol::DESolution,fitness_measure)
    minimum(mapslices(x->fitness_measure(x),sol.u[end],dims = 2))
end

function fitness_evaluation(sol::DESolution,fitness_measure,output::Int64)
    pheno = @view sol.u[end][output,:]
    fitness_measure(pheno)
end

# Evolution 

mutable struct EvolutionaryTrace
    traversed_networks :: Any
    traversed_t2s ::Any
    fitness_trajectory :: Any
    wait_times :: Any
    retcodes :: Any
    converged :: Union{Bool, Vector{Bool}}
    full_weights :: Union{Bool, Vector{Bool}}
    worker_id :: Any
    fitness_transition_times :: Any
    network_transition_times :: Any
    final_networks :: Any
    final_t2s :: Any
    mut_type ::Any
    mut_choices :: Any
    mut_sizes :: Any
end

function has_not_converged(population::Population{Tuple{Float64,Float64}},tolerance::Float64)
    (population.fitness[1] != 0.) || (population.fitness[2] < tolerance)
end

function SSWM_Evolution(start_network::Matrix{Float64},grn_parameters::GRNParameters,β::Union{Float64,Tuple{Float64,Int64}},max_gen::Int64,tolerance::Float64,fitness_function,mutate_function)

    p = (start_network,grn_parameters.degradation)
    
    grn = ODEProblem(gene_regulation_1d!,grn_parameters.g0,(0,Inf),p)

    development = DefaultGRNSolver()
    
    founder = Individual(grn,development)

    founder_fitness = fitness_function(founder.phenotype)

    population = Population(founder,founder_fitness,false)

    gen = 0
    wait_time = 1

    converged = false

    full_weights = false

    evo_trace = EvolutionaryTrace([population.dominant_individual.genotype.p[1]],[population.dominant_individual.phenotype.t[end]],[population.fitness],[],[founder.phenotype.retcode],converged,full_weights,(myid(),gethostname()),[1],[1],[start_network],[population.dominant_individual.phenotype.t[end]],[],[],[])

    while has_not_converged(population,tolerance) && gen < max_gen

        mutant,m_choices,m_type,m_sizes,m_valid = create_mutant(population.dominant_individual,mutate_function,development)

        if m_valid && SciMLBase.successful_retcode(mutant.phenotype.retcode)
            strong_selection!(population,mutant,β,fitness_function)
        else
            population.has_fixed = false
        end

        push!(evo_trace.retcodes,mutant.phenotype.retcode)

        if population.has_fixed
            push!(evo_trace.traversed_networks,population.dominant_individual.genotype.p[1])
            push!(evo_trace.fitness_trajectory,population.fitness)
            push!(evo_trace.wait_times,wait_time)
            push!(evo_trace.traversed_t2s,population.dominant_individual.phenotype.t[end])
            if !isnothing(m_choices)
                push!(evo_trace.mut_choices,m_choices)
                push!(evo_trace.mut_type,m_type)
                push!(evo_trace.mut_sizes,m_sizes)
            end
            wait_time = 1
        else
            wait_time += 1
        end

        gen += 1
    end

    if !has_not_converged(population,tolerance)
        evo_trace.converged = true
    else
        push!(evo_trace.wait_times,wait_time)
    end

    return evo_trace

end

function SSWM_Evolution_SolverIt(start_network::Matrix{Float64},grn_parameters::GRNParameters,solver_it::Int64,β::Union{Float64,Tuple{Float64,Int64}},max_gen::Int64,tolerance::Float64,fitness_function,mutate_function)

    p = (start_network,grn_parameters.degradation)
    
    grn = ODEProblem(gene_regulation_1d!,grn_parameters.g0,(0,Inf),p)

    development = DESystemSolver(Tsit5(),(isoutofdomain=(u,p,t) -> any(x -> x < 0, u), reltol = 1e-6,abstol = 1e-8,callback = TerminateSteadyState(1e-8,1e-6),maxiters = solver_it, verbose = false, save_everystep = false))
    
    founder = Individual(grn,development)

    founder_fitness = fitness_function(founder.phenotype)

    population = Population(founder,founder_fitness,false)

    gen = 0
    wait_time = 1

    converged = false

    full_weights = false

    evo_trace = EvolutionaryTrace([population.dominant_individual.genotype.p[1]],[population.dominant_individual.phenotype.t[end]],[population.fitness],[],[founder.phenotype.retcode],converged,full_weights,(myid(),gethostname()),[1],[1],[start_network],[population.dominant_individual.phenotype.t[end]],[],[],[])

    while has_not_converged(population,tolerance) && gen < max_gen

        mutant,m_choices,m_type,m_sizes,m_valid = create_mutant(population.dominant_individual,mutate_function,development)

        if m_valid && SciMLBase.successful_retcode(mutant.phenotype.retcode)
            strong_selection!(population,mutant,β,fitness_function)
        else
            population.has_fixed = false
        end

        push!(evo_trace.retcodes,mutant.phenotype.retcode)

        if population.has_fixed
            push!(evo_trace.traversed_networks,population.dominant_individual.genotype.p[1])
            push!(evo_trace.fitness_trajectory,population.fitness)
            push!(evo_trace.wait_times,wait_time)
            push!(evo_trace.traversed_t2s,population.dominant_individual.phenotype.t[end])
            if !isnothing(m_choices)
                push!(evo_trace.mut_choices,m_choices)
                push!(evo_trace.mut_type,m_type)
                push!(evo_trace.mut_sizes,m_sizes)
            end
            wait_time = 1
        else
            wait_time += 1
        end

        gen += 1
    end

    if !has_not_converged(population,tolerance)
        evo_trace.converged = true
    else
        push!(evo_trace.wait_times,wait_time)
    end

    return evo_trace

end

# Fitness functions

function malt_fitness_relative(conc,n_stripe::Int64)

    Lt = length(conc)

    N = 2*n_stripe + 1

    id_segments = [Int(floor((k-1)*(Lt/N) + 1)) : Int(floor(k*Lt/N)) for k in 1:N]

    high_sum = 0.
    low_sum = 0.

    for i in 1:N
        if i % 2 == 0
            high_sum += sum(conc[id_segments[i]])
        else
            low_sum += sum(conc[id_segments[i]])
        end
    end

    max_conc = maximum(conc)

    return ((2/(N-1))*high_sum - (2/(N+1))*low_sum) / ((Lt/N)*max_conc)
end

function malt_fitness_absolute(conc,n_stripe::Int64,max_conc::Float64)

    Lt = length(conc)

    N = 2*n_stripe + 1

    id_segments = [Int(floor((k-1)*(Lt/N) + 1)) : Int(floor(k*Lt/N)) for k in 1:N]

    high_sum = 0.
    low_sum = 0.

    for i in 1:N
        if i % 2 == 0
            high_sum += sum(conc[id_segments[i]])
        else
            low_sum += sum(conc[id_segments[i]])
        end
    end

    return ((2/(N-1))*high_sum - (2/(N+1))*low_sum) / ((Lt/N)*max_conc)
end


function malt_fitness_absolute(conc,n_stripe::Int64,max_conc::Float64)

    Lt = length(conc)

    N = 2*n_stripe + 1

    id_segments = [Int(floor((k-1)*(Lt/N) + 1)) : Int(floor(k*Lt/N)) for k in 1:N]

    high_sum = 0.
    low_sum = 0.

    for i in 1:N
        if i % 2 == 0
            high_sum += sum(conc[id_segments[i]])
        else
            low_sum += sum(conc[id_segments[i]])
        end
    end

    return ((2/(N-1))*high_sum - (2/(N+1))*low_sum) / ((Lt/N)*max_conc)
end

function stripe_vary_centre_width(conc,max_conc::Float64,centre_p::Float64,width_p::Float64)

    Lt = length(conc)

    width = width_p*Lt
    centre = centre_p*Lt

    left_boundary = Int(floor(centre - 0.5*width))
    right_boundary = Int(ceil(centre + 0.5*width))

    width_pa = (right_boundary - left_boundary) / Lt

    scale_factor = width/(Lt-width)

    low_sum = scale_factor*(sum(conc[1:left_boundary-1]) + sum(conc[right_boundary+1:Lt]))
    high_sum = sum(conc[left_boundary:right_boundary])


    (1/(width_pa*Lt*max_conc))*(high_sum - low_sum)
end


function stripe_indicator(conc,min_stripe_width,lower_bound,upper_bound)

    low = findall(conc .< lower_bound)
    high = findall(conc .> upper_bound)

    if (length(low) != 0) & (length(high) != 0)

        valid_arrangment = (low[1] < high[1]) & (low[end] > high[end])
        cts_high_region = !(any([(id > high[1]) & (id < high[end]) for id in low]))
        msw = (high[1] - low[1] >= min_stripe_width) & (low[end] - high[end] >= min_stripe_width) & (high[end] - high[1] + 1 >= min_stripe_width)

        if valid_arrangment & cts_high_region & msw
            return 0.
        else
            return -1.
        end
    else
        return -1. 
    end

end

# Dynamical Clustering

function get_rel_dyn_vector(n1,t1,n_steps,save_id)

    w_ind = Individual(n1,t1,grn_parameters,TimeStampedGRNSolver(LinRange(0,t1,n_steps),save_id))

    return reduce(vcat,w_ind.phenotype.u[2:end-1])

end

function get_av_dyn_vector(n1,t1,n_steps,n_segments)

    w_ind = Individual(n1,t1,grn_parameters,TimeStampedGRNSolver(LinRange(0,t1,n_steps)))

    return reduce(vcat,map(x->vec(reduce(hcat,[mean(x[:,t:t+Int(Nc/n_segments)-1],dims = 2) for t in 1:Int(Nc/n_segments):Nc])),w_ind.phenotype.u[2:end-1]))

end

# Minimal Networks

function mask(network,topology)

    new_network = copy(network)

    z0 = findall(x->x == 0,topology)

    new_network[z0] .= 0.

    return new_network

end

function mask_bool(network,mask)

    new_network = copy(network)

    z0 = findall(x->!x,mask)

    new_network[z0] .= 0.

    return new_network

end

function mask_by_id(network,keep_id)

    new_network = copy(network)

    z0 = findall(x->!(x ∈ keep_id),1:length(network))

    new_network[z0] .= 0.

    return new_network

end

function find_minimal_network(network,grn_parameters,development,fitness_function)

    size_S = length(network)

    powerset_S = [[bit == '1' ? true : false for bit in string(i;base = 2,pad = size_S)] for i in 0:2^size_S-1]

    top = []
    top_sizes = []

    ind = Individual(reshape(vcat(network,[0.,0.]),(3,4)),grn_parameters,development)

    for (n,Q) in enumerate(filter(x->x[10] & (x[3] | x[6]),powerset_S))

        size_Q = sum(Q)

        new_network = vcat(mask_bool(network,Q),[0.,0.])

        mutant = Individual(remake(ind.genotype, p = (reshape(new_network,(3,4)),ind.genotype.p[2:end]...)),development)

        mutant_fitness = fitness_function(mutant.phenotype)
        
        if SciMLBase.successful_retcode(mutant.phenotype.retcode) && (abs(mutant_fitness[1]) == 0)
            push!(top,(size_Q,mutant_fitness[2],new_network))
            push!(top_sizes,size_Q)
        end
    end

    if length(top_sizes) > 0

        min_top = minimum(top_sizes)

        return top[findall(top_sizes .== min_top)]

    else
        return top
    end
end

function find_minimal_network_full(network,grn_parameters,development,fitness_function)

    size_S = length(network)

    powerset_S = [[bit == '1' ? true : false for bit in string(i;base = 2,pad = size_S)] for i in 0:2^size_S-1]

    top = []
    top_sizes = []

    ind = Individual(reshape(network,(3,4)),grn_parameters,development)

    for (n,Q) in enumerate(powerset_S)

        size_Q = sum(Q)

        new_network = mask_bool(network,Q)

        mutant = Individual(remake(ind.genotype, p = (reshape(new_network,(3,4)),ind.genotype.p[2:end]...)),development)

        mutant_fitness = fitness_function(mutant.phenotype)
        
        if SciMLBase.successful_retcode(mutant.phenotype.retcode) && (abs(mutant_fitness[1]) == 0)
            push!(top,(size_Q,mutant_fitness[2],new_network))
            push!(top_sizes,size_Q)
        end
    end

    if length(top_sizes) > 0

        min_top = minimum(top_sizes)

        return top[findall(top_sizes .== min_top)]

    else
        return top
    end
end