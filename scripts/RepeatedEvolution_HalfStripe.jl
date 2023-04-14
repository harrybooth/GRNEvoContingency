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

@everywhere include(srcdir("TrajectoryHandling.jl"))
@everywhere include(srcdir("FitnessLandscapes.jl"))

@everywhere start_network_name = "half_right"

# @everywhere example_networks = load(datadir("exp_pro/" * start_network_name * "_networks/examples.jld"))

function repeated_evolution(n_traj,topology = "feed_forward", n_target_stripe =  1, β = Inf, max_gen = 5000, noise_cv = 1., mut_prob = nothing, deletion_prob = 0.05)

    if topology == "feed_forward"
        start_network = [0.0 0.0 0.0 0.28368795845354794; 0.09693796878733349 0.0 0.0 0.0; 0.02660150950444218 -0.26272166357617865 0.6146272196396064 0.0] # right handed feed forward
    elseif topology == "bistable"
        start_network = [0.0 0.0 0.0 0.12728709721871537; -0.014075930837614938 0.0 0.49938778625866675 0.0; 0.0 0.1997636901215515 0.05755668522756788 0.0] # right handed bistable
    elseif topology == "mutual_inh"
        start_network = [0.0 0.0 0.0 -0.06486368943640441; -1.0 0.0 0.7990661288117087 0.0; 0.8329310528122276 0.05802424255402265 0.0 0.0]  # right handed mutual_inh
    elseif topology == "classical"
        start_network = [0.0 0.0 0.0 0.31336321352475677; 0.034733909122486334 0.011312260670141676 0.0 0.0; -0.002116043070043973 -0.1751097063592794 0.5246935289524377 0.0] # right handed classical
    else
        start_network = [0.0 0.0 0.0 0.28368795845354794; 0.09693796878733349 0.0 0.0 0.0; 0.02660150950444218 -0.26272166357617865 0.6146272196396064 0.0] # right handed feed forward
    end

    grn_parameters = DefaultGRNParameters();

    start_ind = Individual(start_network,grn_parameters,DefaultGRNSolver())

    viable_mutations = ones(Int,Ng,Ng+1)
    viable_mutations[2,4] = 0
    viable_mutations[3,4] = 0

    mutation_op = MutationOperator(Normal,(μ = 0.0,σ = noise_cv),findall(viable_mutations .> 0))

    n_sample_func() = rand(Binomial(length(mutation_op.mutation_freq),mut_prob))

    noise_params = (n_sample_func,deletion_prob)
    
    mutate_function = i -> noise(i,mutation_op,noise_params);
    
    output_gene = 3
    
    fitness_function = s -> fitness_evaluation(s,x->malt_fitness(x,n_target_stripe),output_gene);

    # fitness_function = s -> fitness_evaluation(s,x->malt_fitness_left(x),output_gene);

    tolerance = 0.9

    p_final, evo_trace = SSWM_Evolution(start_network,grn_parameters,β,max_gen,tolerance,fitness_function,mutate_function)

    sim = fill((p_final,evo_trace),n_traj)

    @sync for i in 2:n_traj
        @spawn sim[i] = SSWM_Evolution(start_network,grn_parameters,β,max_gen,tolerance,fitness_function,mutate_function)
    end

    print("Evolution complete: " * topology)

    evo_traces_sim = map(x->x[2],sim);

    gt = GenoTrajectories(evo_traces_sim);

    end_fitness = map(x->x[end],gt.fitness_traj);

    converged = end_fitness .> tolerance

    print(string(sum(converged)/length(converged)) * " converged")

    return sim
end

# Make simulation : https://juliadynamics.github.io/DrWatson.jl/dev/workflow/

function makesim(d::Dict)
    
    @unpack n_traj,topology,n_target_stripe, β, max_gen, noise_cv, mut_prob, deletion_prob, start_network_name = d

    sim = repeated_evolution(n_traj,topology,n_target_stripe, β, max_gen, noise_cv, mut_prob, deletion_prob)

    fulld = copy(d)

    # save parameters as part of fulld so that they can be loaded in future notebooks

    fulld["n_max_iters_reached"] = count(x->length(x[2].fitness_trajectory) == max_gen,sim)    
    fulld["describe_proportion_mutants_rejected"] = summarystats(map(x->count(y->y == :MaxIters,x[2].retcodes)/length(x[2].retcodes),sim))

    fulld["raw_data"] = sim

    fulld["parameter_choices"] = Dict("n_traj" => n_traj,"topology" => topologies_test, "n_target_stripe" => n_target_stripe, "β" => β, "max_gen" => max_gen, "noise_cv" => noise_cv, "mut_prob" => mut_prob, "deletion_prob" => deletion_prob, "start_network_name" => start_network_name)

    return fulld
end

# Run

n_traj = 1000
β = 1.
max_gen = 20000
noise_cv = 0.25

n_target_stripe = 1

mut_prob = 0.1
deletion_prob = 0.05

topologies_test = "mutual_inh"

test_specification = Dict("n_traj" => n_traj,"topology" => topologies_test, "n_target_stripe" => n_target_stripe, "β" => β, "max_gen" => max_gen, "noise_cv" => noise_cv, "mut_prob" => mut_prob, "deletion_prob" => deletion_prob,"start_network_name" => start_network_name)

all_tests = dict_list(test_specification);

for (i,d) in enumerate(all_tests)
    f = makesim(d)
    # save(datadir("sims\\repeated_evolution_different_topologies", savename(d, "jld2")), f)
    save(datadir("sims\\repeated_evolution_different_topologies", d["topology"] * "_mi_stripe_1000.jld2"), f)
    # save(datadir("sims/repeated_evolution_different_topologies", "a6.jld2"), f)
end
