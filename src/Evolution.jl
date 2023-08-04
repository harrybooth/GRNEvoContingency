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

struct MutationOperator 
    noise_distribution :: Distribution
    n_sample_func :: Any
    deletion_p :: Float64
    max_w ::Float64
    mutation_weights :: Vector{CartesianIndex{2}}
end

function MutationOperator(noise_distribution,noise_kwargs,n_sample_func,deletion_p,max_w,mutation_freq)
    return MutationOperator(noise_distribution(noise_kwargs...),n_sample_func,deletion_p,max_w,mutation_freq)
end

function create_mutant(ind::Individual,mutate_function,development)
    new_w, m_choices,m_type,m_sizes,valid = mutate_function(ind.genotype.p[1])
    Individual(remake(ind.genotype, p = (new_w,ind.genotype.p[2:end]...)),development),m_choices,m_type,m_sizes,valid
end

function create_mutant(ind::Individual,mutate_function_fm,development,fm_id)
    new_w, m_choices,m_type,m_sizes,valid = mutate_function_fm(ind.genotype.p[1],fm_id)
    Individual(remake(ind.genotype, p = (new_w,ind.genotype.p[2:end]...)),development),m_choices,m_type,m_sizes,valid
end

function create_mutant(ind::Individual,new_w::Matrix{Float64},development)
    Individual(remake(ind.genotype, p = (new_w,ind.genotype.p[2:end]...)),development),nothing,nothing,nothing,nothing
end

function noise(w::Matrix{Float64},mut_op::MutationOperator)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    choices = sample(mut_op.mutation_weights,n_mut,replace = false)

    for index in choices
        if new_w[index] == 0
            # proposal = new_w[index] + rand(mut_op.noise_distribution)
            new_w[index] = rand(Uniform(-mut_op.max_w,mut_op.max_w))
        else
            if rand() < mut_op.deletion_p
                del_index = sample(mut_op.mutation_weights,1)[1]
                new_w[del_index] = 0.

                proposal = new_w[index] + rand(mut_op.noise_distribution)*new_w[index]
                new_w[index] = abs(proposal) > mut_op.max_w ? mut_op.max_w*sign(proposal) : proposal
            else
                proposal = new_w[index] + rand(mut_op.noise_distribution)*new_w[index]
                new_w[index] = abs(proposal) > mut_op.max_w ? mut_op.max_w*sign(proposal) : proposal
            end
        end
    end

    valid = maximum(abs.(new_w)) <= mut_op.max_w

    return new_w,nothing,nothing,nothing, valid
end

function noise_specified(w::Vector{Float64},mut_id::Int64,mut_size::Float64,mut_type::Symbol)
    new_w = copy(w)

    if mut_type == :new
        new_w[mut_id] = mut_size
    elseif mut_type == :existing
        new_w[mut_id] = new_w[mut_id] + mut_size*new_w[mut_id]
    elseif mut_type == :del
        new_w[mut_id] = 0
    else
        new_w[mut_id] = NaN
    end

    return new_w
end

function noise_specified(w::Vector{Float64},mut_id::Vector{Int64},mut_size::Vector{Float64},mut_type::Vector{Any})

    new_w = copy(w)

    for n in 1:length(mut_id)
        if mut_type[n] == :new
            new_w[mut_id[n]] = mut_size[n]
        elseif mut_type[n] == :existing
            new_w[mut_id[n]] = new_w[mut_id[n]] + mut_size[n]*new_w[mut_id[n]]
        elseif mut_type[n] == :del
            new_w[mut_id[n]] = 0
        else
            new_w[mut_id[n]] = NaN
        end
    end

    return new_w
end


function noise_mtype(w::Matrix{Float64},mut_op::MutationOperator)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    choices = sort(sample(mut_op.mutation_weights,n_mut,replace = false))
    mtype = []
    sizes = []

    for index in choices
        if new_w[index] == 0
            n = rand(Uniform(-mut_op.max_w,mut_op.max_w))
            new_w[index] = n
            push!(mtype,:new)
            push!(sizes,n)
        else
            if rand() < mut_op.deletion_p
                push!(sizes,- new_w[index])
                new_w[index] = 0.
                push!(mtype,:del)
            else
                n = rand(mut_op.noise_distribution)
                new_w[index] = new_w[index] + n*new_w[index]
                push!(mtype,:existing)
                push!(sizes,n)
            end
        end
    end

    valid = maximum(abs.(new_w)) <= mut_op.max_w

    return new_w, choices, mtype, sizes, valid
end

function noise_mtype_fm(w::Matrix{Float64},mut_op::MutationOperator,index)

    new_w = copy(w)

    choices = [CartesianIndex(index)]
    mtype = []
    sizes = []

    for index in choices
        if new_w[index] == 0
            n = rand(Uniform(-mut_op.max_w,mut_op.max_w))
            new_w[index] = n
            push!(mtype,:new)
            push!(sizes,n)
        else
            if rand() < mut_op.deletion_p
                push!(sizes,- new_w[index])
                new_w[index] = 0.
                push!(mtype,:del)
            else
                n = rand(mut_op.noise_distribution)
                new_w[index] = new_w[index] + n*new_w[index]
                push!(mtype,:existing)
                push!(sizes,n)
            end
        end
    end

    valid = maximum(abs.(new_w)) <= mut_op.max_w

    return new_w, choices, mtype, sizes, valid
end

function noise_runiform(w::Matrix{Float64},mut_op::MutationOperator,p_uniform)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    choices = sample(mut_op.mutation_weights,n_mut,replace = false)

    for index in choices
        if rand() < p_uniform
            new_w[index] = rand(Uniform(-mut_op.max_w,mut_op.max_w))
        else
            if new_w[index] == 0
                proposal = new_w[index] + rand(mut_op.noise_distribution)
                new_w[index] = abs(proposal) > mut_op.max_w ? mut_op.max_w*sign(proposal) : proposal
            else
                if rand() < mut_op.deletion_p
                    new_w[index] = 0.
                else
                    proposal = new_w[index] + rand(mut_op.noise_distribution)*new_w[index]
                    new_w[index] = abs(proposal) > mut_op.max_w ? mut_op.max_w*sign(proposal) : proposal
                end
            end
        end
    end

    valid = maximum(abs.(new_w)) <= mut_op.max_w

    return new_w,nothing,nothing,nothing, valid
end

function noise_no_additions(w::Matrix{Float64},mut_op::MutationOperator)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    choices = sample(mut_op.mutation_weights,n_mut,replace = false)

    for index in choices

        proposal = new_w[index] + rand(mut_op.noise_distribution)*new_w[index]

        while (sign(proposal) != sign(new_w[index])) || (abs(proposal) > mut_op.max_w)
            proposal = new_w[index] + rand(mut_op.noise_distribution)*new_w[index]
        end

        new_w[index] = proposal
    end

    valid = maximum(abs.(new_w)) <= mut_op.max_w

    return new_w,nothing,nothing,nothing, valid
end

function noise_add(w::Matrix{Float64},mut_op::MutationOperator)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    choices = sample(mut_op.mutation_weights,n_mut,replace = false)

    for index in choices
        proposal = new_w[index] + rand(mut_op.noise_distribution)
        new_w[index] = abs(proposal) > mut_op.max_w ? mut_op.max_w*sign(proposal) : proposal
    end

    valid = maximum(abs.(new_w)) <= mut_op.max_w

    return new_w,nothing,nothing,nothing, valid
end

# Selection 

function fixation_probability(Δf,β)
    1 - exp(-2*β*Δf) 
end

function fixation_probability(Δf1,Δf2,β)
    Δf1 != 0 ? 1 - exp(-2*β*Δf1) : 1 - exp(-2*β*Δf2)
end

function fixation_probability_kim(Δf,β,N)
    (1 - exp(-2*β*Δf)) / (1 - exp(-2*β*N*Δf))
end

function strong_selection!(population::Population{Float64},mutant::Individual,β::Float64,fitness_function)

    mutant_fitness = fitness_function(mutant.phenotype)

    population.has_fixed = false

    if rand() < fixation_probability(mutant_fitness - population.fitness,β)
        population.dominant_individual = mutant
        population.fitness = mutant_fitness
        population.has_fixed = true
    end
end

function strong_selection!(population::Population{Float64},mutant::Individual,β::Tuple{Float64,Int64},fitness_function)

    mutant_fitness = fitness_function(mutant.phenotype)

    has_fixed = false

    if rand() < fixation_probability_kim(mutant_fitness - population.fitness,β[1],β[2])
        population.dominant_individual = mutant
        population.fitness = mutant_fitness
        has_fixed = true
    end
end

function strong_selection!(population::Population{Tuple{Float64,Float64}},mutant::Individual,β::Float64,fitness_function)

    mutant_fitness = fitness_function(mutant.phenotype)

    population.has_fixed = false

    if rand() < fixation_probability(mutant_fitness[1] - population.fitness[1],mutant_fitness[2] - population.fitness[2],β)
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

mutable struct EvoTrace # For legacy data
    traversed_topologies :: Any
    traversed_phenotypes :: Any
    fitness_trajectory :: Any
    retcodes :: Any
end

mutable struct EvolutionaryTrace
    traversed_networks :: Any
    traversed_t2s ::Any
    fitness_trajectory :: Any
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

function has_not_converged(population::Population{Float64},tolerance::Float64)
    population.fitness < tolerance
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

    converged = false

    full_weights = false

    evo_trace = EvolutionaryTrace([population.dominant_individual.genotype.p[1]],[population.dominant_individual.phenotype.t[end]],[population.fitness],[founder.phenotype.retcode],converged,full_weights,(myid(),gethostname()),[1],[1],[start_network],[population.dominant_individual.phenotype.t[end]],[],[],[])

    while has_not_converged(population,tolerance) && gen < max_gen

        mutant,m_choices,m_type,m_sizes,m_valid = create_mutant(population.dominant_individual,mutate_function,development)

        if m_valid && mutant.phenotype.retcode == ReturnCode.Terminated
            strong_selection!(population,mutant,β,fitness_function)
        else
            population.has_fixed = false
        end

        push!(evo_trace.fitness_trajectory,population.fitness)
        push!(evo_trace.retcodes,mutant.phenotype.retcode)

        if population.has_fixed
            push!(evo_trace.traversed_networks,population.dominant_individual.genotype.p[1])
            push!(evo_trace.traversed_t2s,population.dominant_individual.phenotype.t[end])
            if !isnothing(m_choices)
                push!(evo_trace.mut_choices,m_choices)
                push!(evo_trace.mut_type,m_type)
                push!(evo_trace.mut_sizes,m_sizes)
            end
        end

        gen += 1

    end

    if !has_not_converged(population,tolerance)
        evo_trace.converged = true
        final_network = copy(evo_trace.traversed_networks[end])
        if minimum(abs.(final_network[final_network .!= 0.])) > 0.1*maximum(abs.(final_network))  
            evo_trace.full_weights = true
        end
    end

    return evo_trace

end

function SSWM_Evolution_FM(start_network::Matrix{Float64},grn_parameters::GRNParameters,β::Union{Float64,Tuple{Float64,Int64}},max_gen::Int64,tolerance::Float64,fitness_function,mutate_function,mutate_function_fm)

    p = (start_network,grn_parameters.degradation)
    
    grn = ODEProblem(gene_regulation_1d!,grn_parameters.g0,(0,Inf),p)

    development = DefaultGRNSolver()
    
    founder = Individual(grn,development)

    founder_fitness = fitness_function(founder.phenotype)

    population = Population(founder,founder_fitness,false)

    gen = 0

    converged = false

    full_weights = false

    evo_trace = EvolutionaryTrace([population.dominant_individual.genotype.p[1]],[population.dominant_individual.phenotype.t[end]],[population.fitness],[founder.phenotype.retcode],converged,full_weights,(myid(),gethostname()),[1],[1],[start_network],[population.dominant_individual.phenotype.t[end]],[],[],[])

    m_choices_temp = []
    m_type_temp = []
    m_sizes_temp = []

    while (!population.has_fixed) && (gen < max_gen)

        mutant,m_choices,m_type,m_sizes,m_valid = create_mutant(population.dominant_individual,mutate_function_fm,development)

        if m_valid && mutant.phenotype.retcode == ReturnCode.Terminated
            strong_selection!(population,mutant,β,fitness_function)
        else
            population.has_fixed = false
        end

        push!(m_choices_temp,m_choices)
        push!(m_type_temp,m_type)
        push!(m_sizes_temp,m_sizes)

        gen+=1

    end

    push!(evo_trace.fitness_trajectory,population.fitness)
    push!(evo_trace.retcodes,ReturnCode.Terminated)

    push!(evo_trace.traversed_networks,population.dominant_individual.genotype.p[1])
    push!(evo_trace.traversed_t2s,population.dominant_individual.phenotype.t[end])

    if !isnothing(m_choices_temp[end])
        push!(evo_trace.mut_choices,m_choices_temp[end])
        push!(evo_trace.mut_type,m_type_temp[end])
        push!(evo_trace.mut_sizes,m_sizes_temp[end])
    end

    

    while has_not_converged(population,tolerance) && gen < max_gen

        mutant,m_choices,m_type,m_sizes,m_valid = create_mutant(population.dominant_individual,mutate_function,development)

        if m_valid && mutant.phenotype.retcode == ReturnCode.Terminated
            strong_selection!(population,mutant,β,fitness_function)
        else
            population.has_fixed = false
        end

        push!(evo_trace.fitness_trajectory,population.fitness)
        push!(evo_trace.retcodes,mutant.phenotype.retcode)

        if population.has_fixed
            push!(evo_trace.traversed_networks,population.dominant_individual.genotype.p[1])
            push!(evo_trace.traversed_t2s,population.dominant_individual.phenotype.t[end])
            if !isnothing(m_choices)
                push!(evo_trace.mut_choices,m_choices)
                push!(evo_trace.mut_type,m_type)
                push!(evo_trace.mut_sizes,m_sizes)
            end
        end

        gen += 1

    end

    if !has_not_converged(population,tolerance)
        evo_trace.converged = true
        final_network = copy(evo_trace.traversed_networks[end])
        if minimum(abs.(final_network[final_network .!= 0.])) > 0.1*maximum(abs.(final_network))  
            evo_trace.full_weights = true
        end
    end

    return evo_trace

end

function SSWM_Evolution!(population::Population,evo_trace::EvolutionaryTrace,β::Float64,max_gen::Int64,tolerance::Float64,fitness_function,mutate_function)

    development = DefaultGRNSolver()

    gen = 0

    while has_not_converged(population,tolerance) && gen < max_gen

        mutant = create_mutant(population.dominant_individual,mutate_function,development)

        if mutant.phenotype.retcode == ReturnCode.Terminated
            strong_selection!(population,mutant,β,fitness_function)
        else
            population.has_fixed = false
        end

        push!(evo_trace.fitness_trajectory,population.fitness)
        push!(evo_trace.retcodes,mutant.phenotype.retcode)

        if population.has_fixed
            push!(evo_trace.traversed_networks,population.dominant_individual.genotype.p[1])
            push!(evo_trace.traversed_t2s,population.dominant_individual.phenotype.t[end])
        end

        gen += 1

    end

    if !has_not_converged(population,tolerance)
        evo_trace.converged = true
        final_network = copy(evo_trace.traversed_networks[end])
        if minimum(abs.(final_network[final_network .!= 0.])) > 0.1*maximum(abs.(final_network))  
            evo_trace.full_weights = true
        end
    end

end

function SSWM_Evolution!(start_network::Matrix{Float64},evo_trace::EvolutionaryTrace,grn_parameters::GRNParameters,β::Float64,max_gen::Int64,tolerance::Float64,fitness_function,mutate_function)

    p = (start_network,grn_parameters.degradation)
    
    grn = ODEProblem(gene_regulation_1d!,grn_parameters.g0,(0,Inf),p)

    development = DefaultGRNSolver()
    
    founder = Individual(grn,development)

    founder_fitness = fitness_function(founder.phenotype)

    population = Population(founder,founder_fitness,false)

    gen = 0

    evo_trace.converged = false

    while has_not_converged(population,tolerance) && gen < max_gen

        mutant = create_mutant(population.dominant_individual,mutate_function,development)

        if mutant.phenotype.retcode == ReturnCode.Terminated
            strong_selection!(population,mutant,β,fitness_function)
        else
            population.has_fixed = false
        end

        push!(evo_trace.fitness_trajectory,population.fitness)
        push!(evo_trace.retcodes,mutant.phenotype.retcode)

        if population.has_fixed
            push!(evo_trace.traversed_networks,population.dominant_individual.genotype.p[1])
            push!(evo_trace.traversed_t2s,population.dominant_individual.phenotype.t[end])
        end

        gen += 1

    end

    if !has_not_converged(population,tolerance)
        evo_trace.converged = true
        final_network = copy(evo_trace.traversed_networks[end])
        if minimum(abs.(final_network[final_network .!= 0.])) > 0.1*maximum(abs.(final_network))  
            evo_trace.full_weights = true
        end
    end

end

function SSWM_MSelection(start_network::Matrix{Float64},grn_parameters::GRNParameters,β::Union{Float64,Tuple{Float64,Int64}},max_gen::Int64,tolerances::Vector{Float64},fitness_functions::Vector{Function},mutate_function)

    fitness_transition_times = [1]
    network_transition_times = [1]

    final_networks = [start_network]

    converged_status = [true]

    evo_trace = SSWM_Evolution(start_network,grn_parameters,β,max_gen,tolerances[1],fitness_functions[1],mutate_function)

    push!(converged_status,evo_trace.converged)

    push!(fitness_transition_times,length(evo_trace.fitness_trajectory))
    push!(network_transition_times,length(evo_trace.traversed_networks))

    push!(final_networks,evo_trace.traversed_networks[end])

    final_t2s =  []

    push!(final_t2s,evo_trace.traversed_t2s[1])
    push!(final_t2s,evo_trace.traversed_t2s[end])

    for (tol,ff) in zip(tolerances[2:end],fitness_functions[2:end])
        if evo_trace.converged
            SSWM_Evolution!(evo_trace.traversed_networks[end],evo_trace,grn_parameters,β,max_gen,tol,ff,mutate_function)
            push!(converged_status,evo_trace.converged)
            push!(fitness_transition_times,length(evo_trace.fitness_trajectory))
            push!(network_transition_times,length(evo_trace.traversed_networks))
            push!(final_networks,evo_trace.traversed_networks[end])
            push!(final_t2s,evo_trace.traversed_t2s[end])
        else
            push!(converged_status,evo_trace.converged)
            push!(fitness_transition_times,length(evo_trace.fitness_trajectory))
            push!(network_transition_times,length(evo_trace.traversed_networks))
            push!(final_networks,evo_trace.traversed_networks[end])
            push!(final_t2s,evo_trace.traversed_t2s[end])
            break
        end
    end

    evo_trace.converged = converged_status
    evo_trace.fitness_transition_times = fitness_transition_times
    evo_trace.network_transition_times = network_transition_times
    evo_trace.final_networks = final_networks
    evo_trace.final_t2s = final_t2s

    return evo_trace
end

