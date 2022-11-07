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
    DESystemSolver(AutoTsit5(RadauIIA5()),(isoutofdomain=(u,p,t) -> any(x -> x < 0, u), callback = TerminateSteadyState(1e-5,1e-3),maxiters = 1e5, verbose = false, save_everystep = false))
end

# Model parameters

struct GRNParameters
    degradation :: Vector{Float64}
    g0 :: Matrix{Float64}
end

function DefaultGRNParameters()
    GRNParameters(0.05 .* ones(Ng),0.01 .* ones(Ng,Nc))
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

mutable struct Population{T<:Union{Tuple{Float64,Vector{Any}},Tuple{Float64,Float64,Vector{Any}}}}
    dominant_individual::Individual
    fitness :: T
end

function Population(founder::Individual,fitness_function)
    fitness  = fitness_function(founder.phenotype)
    Population(founder,fitness)
end

# Mutation

struct MutationOperator{T} 
    noise_distribution :: Distribution
    mutation_freq :: T
end

function MutationOperator(noise_distribution,noise_kwargs,mutation_freq)
    return MutationOperator(noise_distribution(noise_kwargs...),mutation_freq)
end

function create_mutant(ind::Individual,mutate_function,development)
    Individual(remake(ind.genotype, p = (mutate_function(ind.genotype.p[1]),ind.genotype.p[2:end]...)),development)
end

function create_mutant(ind::Individual,new_w::Matrix{Float64},development)
    Individual(remake(ind.genotype, p = (new_w,ind.genotype.p[2:end]...)),development)
end

function noise(w::Matrix{Float64},mut_op::MutationOperator{Int64},noise_function)
    new_w = copy(w)
    for j in rand(1:size(w,2),mut_op.mutation_freq)
        for i in rand(1:size(w,1),mut_op.mutation_freq)
            new_w[i,j] = noise_function(new_w[i,j],rand(mut_op.noise_distribution))
        end
    end
    return new_w
end

function noise(w::Matrix{Float64},mut_op::MutationOperator{Float64},noise_function)
    new_w = copy(w)
    for j in 1:size(w,2)
        for i in 1:size(w,1)
            if rand() < mut_op.mutation_freq
                new_w[i,j] = noise_function(new_w[i,j],rand(mut_op.noise_distribution))
            end
        end
    end
    return new_w
end

function noise(w::Matrix{Float64},mut_op::MutationOperator{Matrix{Float64}},noise_function)
    new_w = copy(w)
    for j in 1:size(w,2)
        for i in 1:size(w,1)
            if rand() < mut_op.mutation_freq[i,j]
                new_w[i,j] = noise_function(new_w[i,j],rand(mut_op.noise_distribution))
            end
        end
    end
    return new_w
end

function noise(w::Matrix{Float64},mut_op::MutationOperator{Matrix{Int64}},noise_function)
    new_w = copy(w)
    number_of = maximum(mut_op.mutation_freq)
    for index in sample(findall(x-> x > 0,mut_op.mutation_freq),number_of,replace = false)
        new_w[index] = noise_function(new_w[index],rand(mut_op.noise_distribution))
    end
    return new_w
end

# Selection 

function fixation_probability(Δf,β)
    Δf > 0 ? 1 - exp(-2*β*Δf) : 0
end

function fixation_probability(Δf1,Δf2,β)
    if Δf1 > 0 
        1 - exp(-2*β*Δf1)
    elseif (Δf1 == 0) && (Δf2 > 0)
        1 - exp(-2*β*Δf2)
    else
        0
    end
end

function strong_selection!(population::Tuple{Float64,Vector{Any}},mutant::Individual,β::Float64,fitness_function)

    mutant_fitness = fitness_function(mutant.phenotype)

    if rand() < fixation_probability(population.fitness - mutant_fitness,β)
        population.dominant_individual = mutant
        population.fitness = mutant_fitness
    end
end

function strong_selection!(population::Population{Tuple{Float64,Float64,Vector{Any}}},mutant::Individual,β::Float64,fitness_function)

    mutant_fitness = fitness_function(mutant.phenotype)

    if rand() < fixation_probability(population.fitness[1] - mutant_fitness[1],population.fitness[2] - mutant_fitness[2],β)
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

mutable struct EvoTrace
    traversed_topologies :: Any
    fitness_trajectory :: Any
    retcodes :: Any
end

function stopping_criteria(population::Population{Tuple{Float64,Vector{Any}}},tolerance::Float64)
    population.fitness[1] > tolerance
end

function stopping_criteria(population::Population{Tuple{Float64,Float64,Vector{Any}}},tolerance::Float64)
     (population.fitness[2] > tolerance) || (population.fitness[1] != 0)
end

function stopping_criteria(population::Population{Tuple{Float64,Float64,Vector{Any}}},tolerance::Tuple{Float64,Float64})
    (population.fitness[2] > tolerance[2]) || (population.fitness[1] > tolerance[1])
end

function SSWM_Evolution(start_network::Matrix{Float64},grn_parameters::GRNParameters,β::Float64,max_gen::Int64,tolerance::Float64,fitness_function,mutate_function)

    p = (start_network,grn_parameters.degradation)
    
    grn = ODEProblem(gene_regulation_1d!,grn_parameters.g0,(0,Inf),p)

    development = DefaultGRNSolver()
    
    founder = Individual(grn,development)

    population = Population(founder,fitness_function)

    evo_trace = EvoTrace([population.dominant_individual.genotype.p[1]],[population.fitness],[founder.phenotype.retcode])

    gen = 0

    while stopping_criteria(population,tolerance) && gen < max_gen

        mutant = create_mutant(population.dominant_individual,mutate_function,development)

        if mutant.phenotype.retcode == :Terminated
            strong_selection!(population,mutant,β,fitness_function)
        end

        push!(evo_trace.traversed_topologies,population.dominant_individual.genotype.p[1])
        push!(evo_trace.fitness_trajectory,population.fitness)
        push!(evo_trace.retcodes,mutant.phenotype.retcode)

        gen += 1

    end

    return population, evo_trace

end

function SSWM_Evolution_error(start_network::Matrix{Float64},grn_parameters::GRNParameters,β::Float64,max_gen::Int64,tolerance::Float64,fitness_function,mutate_function)

    p = (start_network,grn_parameters.degradation)
    
    grn = ODEProblem(gene_regulation_1d!,grn_parameters.g0,(0,Inf),p)

    development = DefaultGRNSolver()
    
    founder = Individual(grn,development)

    # population = Population(founder,fitness_function)

    # evo_trace = EvoTrace([population.dominant_individual.genotype.p[1]],[population.fitness],[founder.phenotype.retcode])

    # gen = 0

    # w_s = []

    # errors  = []

    # while stopping_criteria(population,tolerance) && gen < max_gen

    #     new_w = mutate_function(population.dominant_individual.genotype.p[1])
    #     push!(w_s,new_w)

    #     try

    #         mutant = create_mutant(population.dominant_individual,new_w,development)

    #         if mutant.phenotype.retcode == :Terminated
    #             strong_selection!(population,mutant,β,fitness_function)
    #         end
    
    #         push!(evo_trace.traversed_topologies,population.dominant_individual.genotype.p[1])
    #         push!(evo_trace.fitness_trajectory,population.fitness)
    #         push!(evo_trace.retcodes,mutant.phenotype.retcode)
    
    #         gen += 1

    #     catch e
    #         push!(errors,gen)
    #         gen += 1
    #     end

    # end

    # return population,w_s,errors

    return founder,fitness_function

end

function SSWM_Evolution!(population::Population,evo_trace::EvoTrace,β::Float64,max_gen::Int64,tolerance::Float64,fitness_function,mutate_function)

    development = DefaultGRNSolver()

    gen = 0

    while stopping_criteria(population,tolerance) && gen < max_gen

        mutant = create_mutant(population.dominant_individual,mutate_function,development)

        if mutant.phenotype.retcode == :Terminated
            strong_selection!(population,mutant,β,fitness_function)
        end

        push!(evo_trace.traversed_topologies,population.dominant_individual.genotype.p[1])
        push!(evo_trace.fitness_trajectory,population.fitness)
        push!(evo_trace.retcodes,mutant.phenotype.retcode)

        gen += 1

    end

end

