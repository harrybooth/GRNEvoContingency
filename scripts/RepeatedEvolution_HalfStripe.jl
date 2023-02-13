using DrWatson

@quickactivate "GRNEvoContingency"


using DifferentialEquations
using Random
using Parameters
using StatsBase
using Printf
using Distributed

using Base.Threads
using Base.Threads: @spawn


@everywhere include(srcdir("Evolution.jl"))
@everywhere include(srcdir("FitnessFunctions.jl"))
@everywhere include(srcdir("TissueModel_ND.jl"))

@everywhere start_network_name = "half_right"

# @everywhere example_networks = load(datadir("exp_pro/" * start_network_name * "_networks/examples.jld"))

function repeated_evolution(n_traj,topology = "classical", n_target_stripe =  1, β = Inf, max_gen = 5000, noise_cv = 1., mut_prob = nothing, deletion_prob = 0.05,noise_method = "additive+deletion", mutation_method = "all_viable")

    start_network = [0.0 0.0 0.0 0.28368795845354794; 0.09693796878733349 0.0 0.0 0.0; 0.02660150950444218 -0.26272166357617865 0.6146272196396064 0.0] # right handed

    grn_parameters = DefaultGRNParameters();

    start_ind = Individual(start_network,grn_parameters,DefaultGRNSolver())

    viable_mutations = ones(Int,Ng,Ng+1)

    mutation_op = MutationOperator(Normal,(μ = 0.0,σ = noise_cv),findall(viable_mutations .> 0))

    n_sample_func() = rand(Binomial(length(mutation_op.mutation_freq),mut_prob))

    noise_params = (n_sample_func,deletion_prob)
    
    mutate_function = i -> noise(i,mutation_op,noise_params);
    
    output_gene = 3
    
    fitness_function = s -> fitness_evaluation(s,x->malt_fitness(x,n_target_stripe),output_gene);

    tolerance = 0.9

    p_final, evo_trace = SSWM_Evolution(start_network,grn_parameters,β,max_gen,tolerance,fitness_function,mutate_function)

    sim = fill((p_final,evo_trace),n_traj)

    @sync for i in 2:n_traj
        @spawn sim[i] = SSWM_Evolution(start_network,grn_parameters,β,max_gen,tolerance,fitness_function,mutate_function)
    end

    return sim
end

# Make simulation : https://juliadynamics.github.io/DrWatson.jl/dev/workflow/

function makesim(d::Dict)
    
    @unpack n_traj,topology,n_target_stripe, β, max_gen, noise_cv, mut_prob, deletion_prob, noise_method, mutation_method,start_network_name = d

    sim = repeated_evolution(n_traj,topology,n_target_stripe, β, max_gen, noise_cv, mut_prob, deletion_prob, noise_method, mutation_method)

    fulld = copy(d)

    # save parameters as part of fulld so that they can be loaded in future notebooks

    fulld["n_max_iters_reached"] = count(x->length(x[2].fitness_trajectory) == max_gen,sim)    
    fulld["describe_proportion_mutants_rejected"] = summarystats(map(x->count(y->y == :MaxIters,x[2].retcodes)/length(x[2].retcodes),sim))

    fulld["raw_data"] = sim

    return fulld
end

# Run

n_traj = 2500
β = Inf
max_gen = 10000
noise_cv = 0.5

n_target_stripe = 1

mut_prob = 0.2
deletion_prob = 0.05

topologies_test = "classical"

test_specification = Dict("n_traj" => n_traj,"topology" => topologies_test, "n_target_stripe" => n_target_stripe, "β" => β, "max_gen" => max_gen, "noise_cv" => noise_cv, "mut_prob" => mut_prob, "deletion_prob" => deletion_prob,
                          "noise_method" => "additive+deletion", "mutation_method" => "all_viable","start_network_name" => start_network_name)

all_tests = dict_list(test_specification);

for (i,d) in enumerate(all_tests)
    f = makesim(d)
    # safesave(datadir("sims/repeated_evolution_different_topologies", savename(d, "jld2")), f)
    save(datadir("sims/repeated_evolution_different_topologies", "a6.jld2"), f)
end
