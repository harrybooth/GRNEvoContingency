using DrWatson

@quickactivate "GRNEvoContingency"

# using Plots

using DifferentialEquations
using Random
using Parameters
using StatsBase
using Printf
using Distributed
using Random
using JLD
using LinearAlgebra
using DataInterpolations
using SparseArrays

using ProgressMeter

@everywhere include(srcdir("TissueModel_ND.jl"))
@everywhere include(srcdir("Evolution.jl"))
@everywhere include(srcdir("FitnessFunctions.jl"))

@everywhere include(srcdir("TrajectoryHandling.jl"))
@everywhere include(srcdir("FitnessLandscapes.jl"))


@everywhere function instability(network_1,network_2,grn_parameters,development,fitness_function)

    itp_ll = InterpolatedLandscapeLean(network_1,network_2,30,grn_parameters,development,fitness_function);

    # next!(progress)

    return 0.5*(itp_ll.slice_fitnesses[1] + itp_ll.slice_fitnesses[end] + 2) - minimum(itp_ll.slice_fitnesses .+ 1)

end

#########################################

grn_parameters = DefaultGRNParameters();

development = DefaultGRNSolver()

output_gene = 3

n_target_stripe = 1

fitness_function = s -> fitness_evaluation(s,x->malt_fitness(x,n_target_stripe),output_gene);

#########################################

evo_traces_cl = load(datadir("sims/repeated_evolution_different_topologies","deletion_prob=0.05_max_gen=40000_mut_prob=0.1_n_target_stripe=1_n_traj=5000_noise_cv=0.5_start_network_name=half_right_topology=classical_β=1.0_1000.jld2"))["data"]
evo_traces_ff = load(datadir("sims/repeated_evolution_different_topologies","deletion_prob=0.05_max_gen=30000_mut_prob=0.1_n_target_stripe=1_n_traj=2000_noise_cv=0.5_start_network_name=half_right_topology=feed_forward_β=1.0_1000.jld2"))["data"]
evo_traces_mi = load(datadir("sims/repeated_evolution_different_topologies","deletion_prob=0.05_max_gen=40000_mut_prob=0.1_n_target_stripe=1_n_traj=5000_noise_cv=0.5_start_network_name=half_right_topology=mutual_inh_β=1.0_1000.jld2"))["data"];

gt_cl = GenoTrajectories(evo_traces_cl);
gt_ff = GenoTrajectories(evo_traces_ff);
gt_mi = GenoTrajectories(evo_traces_mi);

end_networks_cl = map(x->x[:,end],gt_cl.geno_traj);
end_networks_ff = map(x->x[:,end],gt_ff.geno_traj);
end_networks_mi = map(x->x[:,end],gt_mi.geno_traj);

all_networks = reduce(vcat,[end_networks_cl,end_networks_ff,end_networks_mi]);

########################################

n_sample = length(all_networks)

# n_sample = 10

dmat = zeros(n_sample,n_sample)

# p = Progress(Int(0.5*n_sample*n_sample) - n_sample,"Computing pairwise matrix...",50)

@sync for i in 1:n_sample
    for j in 1:n_sample
        if i>j
            @spawn dmat[i,j] = instability(all_networks[i],all_networks[j],grn_parameters,development,fitness_function)
        end
    end
end

# finish!(p)

for i in 1:n_sample
    for j in 1:n_sample
        if i<j
            dmat[i,j] = dmat[j,i]
        end
    end
end

for i in 1:n_sample
    dmat[i,i] = 0.
end

save(datadir("sims/repeated_evolution_different_topologies","pairwise_instabilities_mi_cl_ff_3000.jld2"),"data",dmat)