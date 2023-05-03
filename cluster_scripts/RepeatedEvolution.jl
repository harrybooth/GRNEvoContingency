using DrWatson

cluster_calc = true

if !cluster_calc
    @quickactivate "GRNEvoContingency"
end

using Distributed
using ClusterManagers

# http://jpfairbanks.com/2017/12/27/running-julia-on-slurm-cluster/
# 
# This code is expected to be run from an sbatch script after a module load julia command has been run.
# It starts the remote processes with srun within an allocation specified in the sbatch script.

if cluster_calc
    n_tasks = parse(Int, ENV["SLURM_NTASKS"])
    addprocs(SlurmManager(n_tasks))
    @everywhere using Pkg
    @everywhere Pkg.activate("..")
    @everywhere Pkg.instantiate()
    # @everywhere Pkg.precompile()
end

@everywhere begin
    using DrWatson
    using JLD2
    using Printf
    using Base.Threads
    using Base.Threads: @spawn
end

@everywhere include(srcdir("Evolution.jl"))
@everywhere include(srcdir("FitnessFunctions.jl"))
@everywhere include(srcdir("FitnessLandscapes.jl"))

@everywhere all_experiments = ["RepeatedEvolution_MutualInh"]

@everywhere networks = load(datadir("networks/FindNetworks_HalfStripeLeft_RawData.jld2"));

for exp_name in all_experiments

    @everywhere include(srcdir("ExperimentSetups/" * $exp_name * ".jl"))

    evo_trace = SSWM_Evolution(start_network,grn_parameters,β,max_gen,tolerance,fitness_function,mutate_function)

    evo_trace.converged = false

    sim = fill(evo_trace,n_traj)

    sim_id = findall(x->!x.converged,sim)

    n_trials =  0

    while length(sim_id) != 0

        n_trials += length(sim_id)

        sim = pmap(sim_i -> sim_i.converged ? sim_i : SSWM_Evolution(start_network,grn_parameters,β,max_gen,tolerance,fitness_function,mutate_function),sim)

        sim_id = findall(x->!x.converged,sim)

    end

    summaryd = Dict{String, Any}()

    summaryd["Total traj simulated"] = n_trials
    summaryd["Total traj converged"] = n_traj
    summaryd["ConvergenceRate"] = n_traj/n_trials
    summaryd["N unique workers"] = length(unique(map(et->et.worker_id,sim)))
    summaryd["Average N non-terminated"] = mean(map(et->count(x->x!=:Terminated,et.retcodes) / length(et.retcodes),sim))

    @tag!(summaryd)

    safesave(datadir("exp_summaries",exp_name * "_Summary.jld2"), summaryd)

    ########################################

    geno_traj = map(et->reduce(hcat,map(x->vec(x),et.traversed_networks)),sim);

    @everywhere end_networks = map(x->x[:,end],geno_traj);

    n_sample = length(end_networks)

    print(n_sample)

    dmat = zeros(n_sample,n_sample)

    dmat_id = [(i,j) for i in 1:n_sample, j in 1:n_sample if i>j]

    dmat = pmap(id-> instability_lean(end_networks[id[1]],end_networks[id[2]],N_interp_points,grn_parameters,development,fitness_function),dmat_id)

    ########################################

    fulld = Dict{String, Any}()

    fulld["fitness_traj"] = map(et->et.fitness_trajectory,sim)
    fulld["geno_traj"] = geno_traj
    fulld["retcodes"] = map(et->map(x-> x == :Terminated ? 1 : 0,et.retcodes),sim)
    fulld["dmat"] = dmat
    fulld["dmat_id"] = dmat_id

    @tag!(fulld)

    safesave(datadir("exp_raw",exp_name * "_RawData.jld2"), fulld)

end
