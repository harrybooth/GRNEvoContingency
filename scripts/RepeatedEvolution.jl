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

# Test network function

function repeated_evolution(n_traj, target = [10.,50.,40.], β = Inf, max_gen = 10000, noise_σ = 1., noise_method = "additive", mutation_method = "all_viable")

    start_network =  [ 2.46532  -2.18916   -0.482151  0.100882;
                        0.0       0.336411   0.0       0.0;
                        3.8411   -5.86916    0.0       0.0]


    grn_parameters = DefaultGRNParameters();

    start_ind = Individual(start_network,grn_parameters,DefaultGRNSolver())

    viable_mutations = ones(Ng,Ng+1)
    viable_mutations[2,4] = 0.
    viable_mutations[3,4] = 0.

    mutation_op = MutationOperator(Normal,(μ = 0.0,σ = noise_σ),viable_mutations)

    if noise_method == "multiplicative"
        noise_application = (x,n) -> x * n
    else
        noise_application = (x,n) -> x + n  
    end

    mutate_function = i -> noise(i,mutation_op,noise_application)

    n_stripe = 1

    stripe_threshold = 5.

    min_width = 5.

    output_gene = 3

    fitness_function = s -> fitness_evaluation(s,x->f_sim(x,stripe_threshold,n_stripe,target,min_width),output_gene)

    tolerance = 4.

    p_final, evo_trace = SSWM_Evolution(start_network,grn_parameters,β,max_gen,tolerance,fitness_function,mutate_function)

    sim = fill((p_final,evo_trace),n_traj)

    @sync for i in 2:n_traj
        @spawn sim[i] = SSWM_Evolution(start_network,grn_parameters,β,max_gen,tolerance,fitness_function,mutate_function)
    end

    return sim
end

# Make simulation : https://juliadynamics.github.io/DrWatson.jl/dev/workflow/

function makesim(d::Dict)
    
    @unpack n_traj, target, β, max_gen, noise_σ, noise_method, mutation_method = d

    sim = repeated_evolution(n_traj, target, β, max_gen, noise_σ, noise_method, mutation_method)

    fulld = copy(d)

    fulld["n_max_iters_reached"] = count(x->length(x[2].fitness_trajectory) == max_gen,sim)    
    fulld["describe_proportion_mutants_rejected"] = summarystats(map(x->count(y->y == :MaxIters,x[2].retcodes)/length(x[2].retcodes),sim))

    fulld["raw_data"] = sim

    return fulld
end

# Run

n_traj = 10
β = Inf
max_gen = 10000
noise_σ = 1.

test_specification = Dict("n_traj" => n_traj, "target"=> [target], "β" => β, "max_gen" => max_gen, "noise_σ" => noise_σ, 
                          "noise_method" => "additive", "mutation_method" => "all_viable")

all_tests = dict_list(test_specification);

for (i,d) in enumerate(all_tests)
    f = makesim(d)
    safesave(datadir("sims/repeated_evolution_consistent_start", savename(d, "jld2")), f)
end
