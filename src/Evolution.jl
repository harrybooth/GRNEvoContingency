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

# function Population(founder::Individual,fitness::Float64)
#     # fitness = fitness_function(founder.phenotype)
#     Population(founder,fitness)
# end

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
    Individual(remake(ind.genotype, p = (mutate_function(ind.genotype.p[1]),ind.genotype.p[2:end]...)),development)
end

function create_mutant(ind::Individual,new_w::Matrix{Float64},development)
    Individual(remake(ind.genotype, p = (new_w,ind.genotype.p[2:end]...)),development)
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

    return new_w
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

    return new_w
end

# Selection 

function fixation_probability(Δf,β)
    1 - exp(-2*β*Δf) 
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
    converged :: Bool
    full_weights :: Bool
    worker_id :: Any
end

function has_not_converged(population::Population{Float64},tolerance::Float64)
    population.fitness < tolerance
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

    evo_trace = EvolutionaryTrace([population.dominant_individual.genotype.p[1]],[population.dominant_individual.phenotype.t[end]],[population.fitness],[founder.phenotype.retcode],converged,full_weights,(myid(),gethostname()))

    while has_not_converged(population,tolerance) && gen < max_gen

        mutant = create_mutant(population.dominant_individual,mutate_function,development)

        if mutant.phenotype.retcode == ReturnCode.Terminated
            strong_selection!(population,mutant,β,fitness_function)
        end

        push!(evo_trace.fitness_trajectory,population.fitness)
        push!(evo_trace.retcodes,mutant.phenotype.retcode)

        if population.has_fixed
            push!(evo_trace.traversed_networks,population.dominant_individual.genotype.p[1])
            push!(evo_trace.traversed_t2s,population.dominant_individual.phenotype.t[end])
        end

        gen += 1

    end

    if population.fitness >= tolerance
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
        end

        push!(evo_trace.fitness_trajectory,population.fitness)
        push!(evo_trace.retcodes,mutant.phenotype.retcode)

        if has_fixed
            push!(evo_trace.traversed_networks,population.dominant_individual.genotype.p[1])
            push!(evo_trace.traversed_t2s,population.dominant_individual.phenotype.t[end])
        end

        gen += 1

    end

    if population.fitness >= tolerance
        evo_trace.converged = true
    end

end

