using DifferentialEquations
using Distributions
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

# Model parameters

struct GRNParameters
    degradation :: Vector{Float64}
    g0 :: Matrix{Float64}
end

function DefaultGRNParameters()
    GRNParameters(degradation_rate .* ones(Ng),initial_concentration .* ones(Ng,Nc))
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

mutable struct Population{T}
    dominant_individual::Individual
    fitness :: T
end

function Population(founder::Individual,fitness_function)
    fitness = fitness_function(founder.phenotype)
    Population(founder,fitness)
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

# Selection 

function fixation_probability(Δf,β)
    1 - exp(-2*β*Δf) 
end

function fixation_probability_kim(Δf,β,N)
    (1 - exp(-2*β*Δf)) / (1 - exp(-2*β*N*Δf))
end

function strong_selection!(population::Population{Float64},mutant::Individual,β::Float64,fitness_function)

    mutant_fitness = fitness_function(mutant.phenotype)

    if rand() < fixation_probability(mutant_fitness - population.fitness,β)
        population.dominant_individual = mutant
        population.fitness = mutant_fitness
    end
end

function strong_selection!(population::Population{Float64},mutant::Individual,β::Tuple{Float64,Int64},fitness_function)

    mutant_fitness = fitness_function(mutant.phenotype)

    if rand() < fixation_probability_kim(mutant_fitness - population.fitness,β[1],β[2])
        population.dominant_individual = mutant
        population.fitness = mutant_fitness
    end
end

# Fitness fitness_evaluation

function fitness_evaluation(sol::DESolution,fitness_measure)
    minimum(mapslices(x->fitness_measure(x),sol.u[end],dims = 2))
end

function fitness_evaluation(sol::DESolution,fitness_measure,output::Int64)
    fitness_measure(sol.u[end][output,:])
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
    fitness_trajectory :: Any
    retcodes :: Any
end

function stopping_criteria(population::Population{Float64},tolerance::Float64)
    population.fitness < tolerance
end

function SSWM_Evolution(start_network::Matrix{Float64},grn_parameters::GRNParameters,β::Union{Float64,Tuple{Float64,Int64}},max_gen::Int64,tolerance::Float64,fitness_function,mutate_function)

    p = (start_network,grn_parameters.degradation)
    
    grn = ODEProblem(gene_regulation_1d!,grn_parameters.g0,(0,Inf),p)

    development = DefaultGRNSolver()
    
    founder = Individual(grn,development)

    population = Population(founder,fitness_function)

    evo_trace = EvolutionaryTrace([population.dominant_individual.genotype.p[1]],[population.fitness],[founder.phenotype.retcode])

    gen = 0

    while stopping_criteria(population,tolerance) && gen < max_gen

        mutant = create_mutant(population.dominant_individual,mutate_function,development)

        if mutant.phenotype.retcode == :Terminated
            strong_selection!(population,mutant,β,fitness_function)
        end

        push!(evo_trace.fitness_trajectory,population.fitness)
        push!(evo_trace.retcodes,mutant.phenotype.retcode)
        push!(evo_trace.traversed_networks,population.dominant_individual.genotype.p[1])

        gen += 1

    end

    return evo_trace

end

function SSWM_Evolution!(population::Population,evo_trace::EvolutionaryTrace,β::Float64,max_gen::Int64,tolerance::Float64,fitness_function,mutate_function)

    development = DefaultGRNSolver()

    gen = 0

    while stopping_criteria(population,tolerance) && gen < max_gen

        mutant = create_mutant(population.dominant_individual,mutate_function,development)

        if mutant.phenotype.retcode == :Terminated
            strong_selection!(population,mutant,β,fitness_function)
        end

        push!(evo_trace.fitness_trajectory,population.fitness)
        push!(evo_trace.retcodes,mutant.phenotype.retcode)
        push!(evo_trace.traversed_topologies,population.dominant_individual.genotype.p[1])

        gen += 1

    end

end

