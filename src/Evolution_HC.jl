using HomotopyContinuation
using Symbolics
using LinearAlgebra
using Distributions
using StatsBase
using Random

# Model

@var w[1:3, 1:4]  g[1:3], d[1:3], m

s(i) = sum(w[i,k] * g[k] for k in 1:3) + w[i, 4] * m

function Δ(i)
  I = s(i) + h_a
  [
    # d[j,i] = sqrt(I^2 + 1)
    d[i]^2 - (I^2 + h_b)
    # σ(s(i,j)) - λ_g g_j^{(i)}
    0.5 * ((I / d[i]) + 1) - λ_g * g[i]
  ]
end

function GRN()
  System(
      [Δ(1)
       Δ(2)
       Δ(3)];
    variables=[vec(g);vec(d)],
    parameters=[vec(w);m])
end

# Symbolics - create jacobian

@variables G[1:3],M

@variables W[1:3,1:4]

g_sum(i) = sum(W[i,k] * G[k] for k in 1:3) + W[i, 4] * M

g_delta(i) = σ(g_sum(i)) - λ_g*G[i]

J = Symbolics.jacobian([g_delta(1),g_delta(2),g_delta(3)], G)

J_call = eval(Symbolics.build_function(J,[G;vec(W); M])[1])

# Define functions which will specify whether a solution is valid or not (positive real +  all real eigenvalues negative)

function satisfies_system(G,p)

    G1 = G[1]
    G2 = G[2]
    G3 = G[3]

    w11, w21, w31, w12, w22, w32, w13, w23, w33, w1m, w2m, w3m, M = p

    r1 = σ(w11*G1 + w12*G2 + w13*G3 + w1m*M) - λ_g*G1
    r2 = σ(w21*G1 + w22*G2 + w23*G3 + w2m*M) - λ_g*G2
    r3 = σ(w31*G1 + w32*G2 + w33*G3 + w3m*M) - λ_g*G3

    eig = J_call([G1,G2,G3,w11, w21, w31, w12, w22, w32, w13, w23, w33, w1m, w2m, w3m, M])

    return all(abs.([r1,r2,r3]) .< 1e-11) && all(real.(eigvals(eig)) .< 0)
end

function get_valid_solutions(r,p)
    rp = filter(s -> all(s .> 0), real_solutions(r))
    rpv = filter(x->satisfies_system(x,p),rp)
    return rpv
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
    system :: System
    genotype :: Matrix{Float64}
    phenotype :: Matrix{Float64}
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
    pheno_class :: Any
end

function Population(founder::Individual,fitness_function)
    fitness,pheno_class  = fitness_function(founder.phenotype)
    Population(founder,fitness,pheno_class)
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

# function noise(w::Matrix{Float64},mut_op::MutationOperator{Matrix{Float64}},noise_function)
#     new_w = copy(w)
#     for j in 1:size(w,2)
#         for i in 1:size(w,1)
#             if rand() < mut_op.mutation_freq[i,j]
#                 new_w[i,j] = noise_function(new_w[i,j],rand(mut_op.noise_distribution))
#             end
#         end
#     end
#     return new_w
# end

function noise(w::Matrix{Float64},mut_op::MutationOperator{Matrix{Int64}},noise_function)
    new_w = copy(w)
    # number_of = maximum(mut_op.mutation_freq)
    for index in sample(findall(x-> x > 0,mut_op.mutation_freq),1,replace = false)
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

function strong_selection!(population::Population{Float64},mutant::Individual,β::Float64,fitness_function)

    mutant_fitness,mutant_pheno_class = fitness_function(mutant.phenotype)

    if rand() < fixation_probability(mutant_fitness - population.fitness,β)
        population.dominant_individual = mutant
        population.fitness = mutant_fitness
        population.pheno_class = mutant_pheno_class
    end
end

function strong_selection!(population::Population{Tuple{Float64,Float64}},mutant::Individual,β::Float64,fitness_function)

    mutant_fitness,mutant_pheno_class = fitness_function(mutant.phenotype)

    if rand() < fixation_probability(mutant_fitness[1] - population.fitness[1],mutant_fitness[2] - population.fitness[2],β)
        population.dominant_individual = mutant
        population.fitness = mutant_fitness
        population.pheno_class = mutant_pheno_class
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
    traversed_phenotypes :: Any
    fitness_trajectory :: Any
    retcodes :: Any
end

function stopping_criteria(population::Population{Float64},tolerance::Float64)
    population.fitness > tolerance
end

function stopping_criteria(population::Population{Tuple{Float64,Float64}},tolerance::Float64)
     (population.fitness[2] > tolerance) || (population.fitness[1] != 0)
end

function stopping_criteria(population::Population{Tuple{Float64,Float64}},tolerance::Tuple{Float64,Float64})
    (population.fitness[2] > tolerance[2]) || (population.fitness[1] > tolerance[1])
end

function SSWM_Evolution(start_network::Matrix{Float64},grn_parameters::GRNParameters,β::Float64,max_gen::Int64,tolerance::Float64,fitness_function,mutate_function)

    p = (start_network,grn_parameters.degradation)
    
    grn = ODEProblem(gene_regulation_1d!,grn_parameters.g0,(0,Inf),p)

    development = DefaultGRNSolver()
    
    founder = Individual(grn,development)

    population = Population(founder,fitness_function)

    evo_trace = EvoTrace([population.dominant_individual.genotype.p[1]],[population.pheno_class],[population.fitness],[founder.phenotype.retcode])

    gen = 0

    while stopping_criteria(population,tolerance) && gen < max_gen

        mutant = create_mutant(population.dominant_individual,mutate_function,development)

        if mutant.phenotype.retcode == :Terminated
            strong_selection!(population,mutant,β,fitness_function)
        end

        push!(evo_trace.traversed_topologies,population.dominant_individual.genotype.p[1])
        push!(evo_trace.traversed_phenotypes,population.pheno_class)
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
        push!(evo_trace.traversed_phenotypes,population.pheno_class)
        push!(evo_trace.fitness_trajectory,population.fitness)
        push!(evo_trace.retcodes,mutant.phenotype.retcode)

        gen += 1

    end

end

