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
@everywhere include(srcdir("NetworkTopologies.jl"))

# Test network function


function new_start(grn_parameters,β,max_gen,tolerance,fitness_function,noise_σ)

    choice = rand(1:6)

    start_topology = network_topology_list[choice]

    mutation_op = MutationOperator(Normal,(μ = 0.0,σ = noise_σ),start_topology)

    noise_application = (x,n) -> x + n  

    mutate_function = i -> noise(i,mutation_op,noise_application)

    p_final, evo_trace = SSWM_Evolution((0.9995 .^ rand(0:10000,Ng,Ng+1)) .* 10 .* rand(Ng,Ng+1) .* start_topology,grn_parameters,β,max_gen,tolerance,fitness_function,mutate_function)

    return choice, p_final, evo_trace
end


function find_topologies(n_traj, target = [85.,30.,85.], β = Inf, max_gen = 10000, noise_σ = 1.)

    grn_parameters = DefaultGRNParameters()

    n_stripe = 1

    stripe_threshold = 5.

    min_width = 5.

    output_gene = 3

    fitness_function = s -> fitness_evaluation(s,x->f_sim(x,stripe_threshold,n_stripe,target,min_width),output_gene)

    tolerance = 4.

    choice, p_final, evo_trace = new_start(grn_parameters,β,max_gen,tolerance,fitness_function,noise_σ)

    sim = fill((choice,p_final,evo_trace),n_traj)

    @sync for i in 2:n_traj
        @spawn sim[i] = new_start(grn_parameters,β,max_gen,tolerance,fitness_function,noise_σ)
    end

    return sim
end

# Make simulation : https://juliadynamics.github.io/DrWatson.jl/dev/workflow/

function makesim(d::Dict)
    
    @unpack n_traj, target, β, max_gen, noise_σ  = d

    sim = find_topologies(n_traj, target, β, max_gen, noise_σ)

    fulld = copy(d)

    fulld["n_max_gen_reached"] = count(x->length(x[3].fitness_trajectory) == max_gen,sim)    
    fulld["describe_proportion_mutants_rejected"] = summarystats(map(x->count(y->y == :MaxIters,x[3].retcodes)/length(x[3].retcodes),sim))

    fulld["best_fitness"] = minimum(map(x->x[3].fitness_trajectory[2],sim))

    fulld["raw_data"] = sim

    return fulld
end

# Run

n_traj = 1
β = Inf
max_gen = 10000
noise_σ = 1.

target = [80.,40.,80.]

test_specification = Dict("n_traj" => n_traj, "target"=> [target], "β" => β, "max_gen" => max_gen, "noise_σ" => noise_σ)

all_tests = dict_list(test_specification);

for (i,d) in enumerate(all_tests)
    f = makesim(d)
    safesave(datadir("sims/find_topologies", savename(d, "jld2")), f)
end
