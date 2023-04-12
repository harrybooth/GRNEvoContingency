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

using JLD
using LinearAlgebra
using DataInterpolations

@everywhere include(srcdir("Evolution.jl"))
@everywhere include(srcdir("FitnessFunctions.jl"))
@everywhere include(srcdir("TissueModel_ND.jl"))

@everywhere include(srcdir("TrajectoryHandling.jl"))
@everywhere include(srcdir("FitnessLandscapes.jl"))


@everywhere function instability(network_1,network_2,grn_parameters,development,fitness_function)

    itp_ll = InterpolatedLandscapeLean(network_1,network_2,30,grn_parameters,development,fitness_function);

    return 0.5*(itp_ll.slice_fitnesses[1] + itp_ll.slice_fitnesses[end] + 2) - minimum(itp_ll.slice_fitnesses .+ 1)

end

# @everywhere example_networks = load(datadir("exp_pro/" * start_network_name * "_networks/examples.jld"))

function repeated_evolution(all_networks, all_labels, β = Inf, max_gen = 5000, noise_cv = 1., mut_prob = nothing, deletion_prob = 0.05)

    n_traj = length(all_networks)

    grn_parameters = DefaultGRNParameters();

    development = DefaultGRNSolver()

    start_ind = Individual(reshape(all_networks[1],(3,4)),grn_parameters,development)

    viable_mutations = ones(Int,Ng,Ng+1)

    mutation_op = MutationOperator(Normal,(μ = 0.0,σ = noise_cv),findall(viable_mutations .> 0))

    n_sample_func() = rand(Binomial(length(mutation_op.mutation_freq),mut_prob))

    noise_params = (n_sample_func,deletion_prob)
    
    mutate_function = i -> noise(i,mutation_op,noise_params);
    
    output_gene = 3

    fitness_function_start = s -> fitness_evaluation(s,x->malt_fitness_left(x),output_gene);
    
    fitness_function = s -> fitness_evaluation(s,x->malt_fitness(x,1),output_gene);

    tolerance = 0.9

    p_final, evo_trace = SSWM_Evolution(reshape(all_networks[1],(3,4)),grn_parameters,β,max_gen,tolerance,fitness_function,mutate_function)

    sim = fill((p_final,evo_trace),n_traj)

    @sync for i in 2:n_traj
        @spawn sim[i] = SSWM_Evolution(reshape(all_networks[i],(3,4)),grn_parameters,β,max_gen,tolerance,fitness_function,mutate_function)
    end

    print("Evolution complete: ")

    evo_traces_sim = map(x->x[2],sim);

    gt = GenoTrajectories(evo_traces_sim);

    end_fitness = map(x->x[end],gt.fitness_traj);

    converged = end_fitness .> tolerance

    print("Evolution complete: " * string(sum(converged)) * " converged")

    start_networks = map(x->x[:,1],gt.geno_traj[converged]);

    end_networks = map(x->x[:,end],gt.geno_traj[converged]);

    final_labels = all_labels[converged]

    n_sample = length(end_networks)

    dmat_start = zeros(n_sample,n_sample)

    dmat_id = [(i,j) for i in 1:n_sample, j in 1:n_sample if i>j]

    @sync for id in dmat_id
        @spawn dmat_start[id...] = instability(start_networks[id[1]],start_networks[id[2]],grn_parameters,development,fitness_function_start)
    end

    print("Pairwise instability complete (start): ")

    dmat_end = zeros(n_sample,n_sample)

    @sync for id in dmat_id
        @spawn dmat_end[id...] = instability(end_networks[id[1]],end_networks[id[2]],grn_parameters,development,fitness_function)
    end

    print("Pairwise instability complete (end): ")

    return sim[converged],dmat_start, dmat_end, final_labels
end

# Make simulation : https://juliadynamics.github.io/DrWatson.jl/dev/workflow/

function makesim(d::Dict)

    # evo_traces_cl = load(datadir("sims/repeated_evolution_different_topologies","deletion_prob=0.05_max_gen=40000_mut_prob=0.1_n_target_stripe=1_n_traj=5000_noise_cv=0.5_start_network_name=half_right_topology=classical_β=1.0_1000.jld2"))["data"]
    # evo_traces_ff = load(datadir("sims/repeated_evolution_different_topologies","deletion_prob=0.05_max_gen=30000_mut_prob=0.1_n_target_stripe=1_n_traj=2000_noise_cv=0.5_start_network_name=half_right_topology=feed_forward_β=1.0_1000.jld2"))["data"]
    # evo_traces_mi = load(datadir("sims/repeated_evolution_different_topologies","deletion_prob=0.05_max_gen=40000_mut_prob=0.1_n_target_stripe=1_n_traj=5000_noise_cv=0.5_start_network_name=half_right_topology=mutual_inh_β=1.0_1000.jld2"))["data"];

    # n_sub_sample = 500

    # gt_cl = GenoTrajectories(evo_traces_cl[1:n_sub_sample]);
    # gt_ff = GenoTrajectories(evo_traces_ff[1:n_sub_sample]);
    # gt_mi = GenoTrajectories(evo_traces_mi[1:n_sub_sample]);

    run_data_cl = load(datadir("sims/repeated_evolution_different_topologies","classical_halfright_to_halfleft.jld2"))["raw_data"]
    run_data_ff = load(datadir("sims/repeated_evolution_different_topologies","bistable_halfright_to_halfleft.jld2"))["raw_data"]
    run_data_mi = load(datadir("sims/repeated_evolution_different_topologies","mutual_inh_halfright_to_halfleft.jld2"))["raw_data"];

    evo_traces_cl = map(x->x[2],run_data_cl);
    evo_traces_ff = map(x->x[2],run_data_ff);
    evo_traces_mi = map(x->x[2],run_data_mi);

    # n_sub_sample = 5

    # gt_cl = GenoTrajectories(evo_traces_cl[1:n_sub_sample]);
    # gt_ff = GenoTrajectories(evo_traces_ff[1:n_sub_sample]);
    # gt_mi = GenoTrajectories(evo_traces_mi[1:n_sub_sample]);

    gt_cl = GenoTrajectories(evo_traces_cl);
    gt_ff = GenoTrajectories(evo_traces_ff);
    gt_mi = GenoTrajectories(evo_traces_mi);

    end_fitness_cl = map(x->x[end],gt_cl.fitness_traj);
    end_fitness_ff = map(x->x[end],gt_ff.fitness_traj);
    end_fitness_mi = map(x->x[end],gt_mi.fitness_traj);

    converged_cl = end_fitness_cl .> 0.9
    converged_ff = end_fitness_ff .> 0.9
    converged_mi = end_fitness_mi .> 0.9

    end_networks_cl = map(x->x[:,end],gt_cl.geno_traj[converged_cl]);
    end_networks_ff = map(x->x[:,end],gt_ff.geno_traj[converged_ff]);
    end_networks_mi = map(x->x[:,end],gt_mi.geno_traj[converged_mi]);

    label_cl = [1 for n in end_networks_cl]
    label_ff = [2 for n in end_networks_ff]
    label_mi = [3 for n in end_networks_mi]

    all_labels = reduce(vcat,[label_cl,label_ff,label_mi]);

    all_networks = reduce(vcat,[end_networks_cl,end_networks_ff,end_networks_mi]);

    print(length(all_networks))
    
    @unpack β, max_gen, noise_cv, mut_prob, deletion_prob = d

    sim,dmat_start, dmat_end, final_labels = repeated_evolution(all_networks,all_labels, β, max_gen, noise_cv, mut_prob, deletion_prob)

    fulld = Dict()

    # save parameters as part of fulld so that they can be loaded in future notebooks

    fulld["raw_data"] = sim
    fulld["dmat_start"] = dmat_start
    fulld["dmat_end"] = dmat_end
    fulld["labels"] = final_labels
    fulld["parameters"] = d

    return fulld
end

# Run

β = 1.
max_gen = 40000
noise_cv = 0.25

mut_prob = 0.1
deletion_prob = 0.05

topologies_test = ["mutual_inh","classical"]

test_specification = Dict("β" => β, "max_gen" => max_gen, "noise_cv" => noise_cv, "mut_prob" => mut_prob, "deletion_prob" => deletion_prob)

all_tests = dict_list(test_specification);

for (i,d) in enumerate(all_tests)
    f = makesim(d)
    # save(datadir("sims\\repeated_evolution_different_topologies", savename(d, "jld2")), f)
    save(datadir("sims/repeated_evolution_different_topologies", "aim_stripe_from_left.jld2"), f)
end
