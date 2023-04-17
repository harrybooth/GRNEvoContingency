using DrWatson

@quickactivate "GRNEvoContingency"

# using Plots

# using Plots
using Random
using Printf
using JLD
using LinearAlgebra
using NearestNeighbors
using DataInterpolations
using Clustering
using ClusterValidityIndices
using MultivariateStats
using KernelDensity
using NetworkLayout
using TravelingSalesmanExact, GLPK

using Distributed

using Memoization
using SparseArrays

using Base.Threads
using Base.Threads: @spawn

@everywhere include(srcdir("TissueModel_ND.jl"))
@everywhere include(srcdir("Evolution.jl"))
@everywhere include(srcdir("FitnessFunctions.jl"))

@everywhere include(srcdir("TrajectoryHandling.jl"))
@everywhere include(srcdir("FitnessLandscapes.jl"))


@everywhere function instability(network_1,network_2,grn_parameters,development,fitness_function)

    itp_ll = InterpolatedLandscapeLean(network_1,network_2,30,grn_parameters,development,fitness_function);

    return 0.5*(itp_ll.slice_fitnesses[1] + itp_ll.slice_fitnesses[end] + 2) - minimum(itp_ll.slice_fitnesses .+ 1)

end

#########################################

@everywhere grn_parameters = DefaultGRNParameters();

@everywhere development = DefaultGRNSolver()

@everywhere output_gene = 3

@everywhere n_target_stripe = 1

@everywhere fitness_function = s -> fitness_evaluation(s,x->malt_fitness(x,n_target_stripe),output_gene);

#########################################

# evo_traces_cl = load(datadir("sims/repeated_evolution_different_topologies","deletion_prob=0.05_max_gen=40000_mut_prob=0.1_n_target_stripe=1_n_traj=5000_noise_cv=0.5_start_network_name=half_right_topology=classical_β=1.0_1000.jld2"))["data"]
# evo_traces_ff = load(datadir("sims/repeated_evolution_different_topologies","deletion_prob=0.05_max_gen=30000_mut_prob=0.1_n_target_stripe=1_n_traj=2000_noise_cv=0.5_start_network_name=half_right_topology=feed_forward_β=1.0_1000.jld2"))["data"]
# evo_traces_mi = load(datadir("sims/repeated_evolution_different_topologies","mutual_inhmi_stripe_0.5_mutation_rate_1000.jld2"))["raw_data"];

run_data_mi = load(datadir("sims/repeated_evolution_different_topologies","mutual_inhmi_stripe_0.5_mutation_rate_1000.jld2"))["raw_data"];
evo_traces_mi = map(x->x[2],run_data_mi);


gt_mi = GenoTrajectories(evo_traces_mi);

end_fitness_mi = map(x->x[end],gt_mi.fitness_traj);

converged_mi = end_fitness_mi .> 0.9

end_networks_mi = map(x->x[:,end],gt_mi.geno_traj[converged_mi]);

n_sub_sample = 1000

# gt_cl = GenoTrajectories(evo_traces_cl[1:n_sub_sample]);
# gt_ff = GenoTrajectories(evo_traces_ff[1:n_sub_sample]);
# gt_mi = GenoTrajectories(evo_traces_mi[1:n_sub_sample]);

# end_networks_cl = map(x->x[:,end],gt_cl.geno_traj);
# end_networks_ff = map(x->x[:,end],gt_ff.geno_traj);
# end_networks_mi = map(x->x[:,end],gt_mi.geno_traj);

# all_networks = reduce(vcat,[end_networks_cl,end_networks_ff,end_networks_mi]);

all_networks = copy(end_networks_mi[1:n_sub_sample])

print("loaded data...")

########################################

n_sample = length(all_networks)

print(n_sample)

dmat = zeros(n_sample,n_sample)

dmat_id = [(i,j) for i in 1:n_sample, j in 1:n_sample if i>j]

@sync for id in dmat_id
    @spawn dmat[id...] = instability(all_networks[id[1]],all_networks[id[2]],grn_parameters,development,fitness_function)
end

print("finished comp...") 

# finish!(p)

# for i in 1:n_sample
#     for j in 1:n_sample
#         if i<j
#             dmat[i,j] = dmat[j,i]
#         end
#     end
# end

# for i in 1:n_sample
#     dmat[i,i] = 0.
# end

# save(datadir("sims/repeated_evolution_different_topologies","pairwise_instabilities_mi_1000_" * string(n_sample) * ".jld2"),"data",dmat)

save(datadir("sims/repeated_evolution_different_topologies","mutual_inhmi_stripe_0.5_mutation_rate_1000_pairwise_lmc.jld2"),"data",dmat)