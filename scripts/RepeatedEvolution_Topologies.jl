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

@everywhere example_networks = load(datadir("exp_pro/80-40-80_networks/examples.jld"))

function repeated_evolution(topology,n_traj, target =  [(40.,20.)], β = Inf, max_gen = 20000, noise_cv = 1., noise_method = "additive", mutation_method = "all_viable")

    start_network = example_networks[topology]

    grn_parameters = DefaultGRNParameters();

    start_ind = Individual(start_network,grn_parameters,DefaultGRNSolver())

    viable_mutations = ones(Int,Ng,Ng+1)
    viable_mutations[2,4] = 0
    viable_mutations[3,4] = 0

    mutation_op = MutationOperator(Normal,(μ = 0.0,σ = noise_cv),viable_mutations)

    if noise_method == "multiplicative"
        noise_application = (x,n) -> x * n
    else
        noise_application = (x,n) -> x + x*n  
    end

    mutate_function = i -> noise(i,mutation_op,noise_application)

    n_stripe = 1

    stripe_threshold = 5.

    min_width = 5.

    output_gene = 3

    fitness_function = s -> fitness_evaluation(s,x->f_sim_cw(x,stripe_threshold,n_stripe,target,min_width),output_gene)

    tolerance = 2.

    p_final, evo_trace = SSWM_Evolution(start_network,grn_parameters,β,max_gen,tolerance,fitness_function,mutate_function)

    sim = fill((p_final,evo_trace),n_traj)

    @sync for i in 2:n_traj
        @spawn sim[i] = SSWM_Evolution(start_network,grn_parameters,β,max_gen,tolerance,fitness_function,mutate_function)
    end

    return sim
end

# Make simulation : https://juliadynamics.github.io/DrWatson.jl/dev/workflow/

function makesim(d::Dict)
    
    @unpack topology,n_traj, target, β, max_gen, noise_cv, noise_method, mutation_method = d

    sim = repeated_evolution(topology,n_traj, target, β, max_gen, noise_cv, noise_method, mutation_method)

    fulld = copy(d)

    fulld["n_max_iters_reached"] = count(x->length(x[2].fitness_trajectory) == max_gen,sim)    
    fulld["describe_proportion_mutants_rejected"] = summarystats(map(x->count(y->y == :MaxIters,x[2].retcodes)/length(x[2].retcodes),sim))

    fulld["raw_data"] = sim

    return fulld
end

# Run

n_traj = 100
β = Inf
max_gen = 20000
noise_cv = 1.

target = [[(40.,20.)]]

topologies_test = collect(keys(filter(x->typeof(x[2]) == Matrix{Float64},example_networks)))

test_specification = Dict("topology" => topologies_test, "n_traj" => n_traj, "target"=> target, "β" => β, "max_gen" => max_gen, "noise_cv" => noise_cv, 
                          "noise_method" => "additive", "mutation_method" => "all_viable")

all_tests = dict_list(test_specification);

for (i,d) in enumerate(all_tests)
    f = makesim(d)
    safesave(datadir("sims/repeated_evolution_different_topologies", savename(d, "jld2")), f)
end
