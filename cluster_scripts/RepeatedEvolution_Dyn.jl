# This code is expected to be run from an sbatch script after a module load julia command has been run.
# It starts the remote processes with srun within an allocation specified in the sbatch script.

using Pkg

Pkg.activate("..")
Pkg.instantiate()
Pkg.precompile()

using DrWatson
using Distributed
using ClusterManagers

projectdir_static = dirname(Base.active_project())

cluster_calc = true

if !cluster_calc
    @quickactivate "GRNEvoContingency"
end

if cluster_calc
    n_tasks = parse(Int, ENV["SLURM_NTASKS"])
    addprocs(SlurmManager(n_tasks))
    @everywhere using Pkg
    @everywhere Pkg.activate("..")
    @everywhere Pkg.instantiate()
    @everywhere Pkg.precompile()
end

@everywhere begin
    using DrWatson
    using JLD2
    using Printf
    using Base.Threads
    using Base.Threads: @spawn
end

@everywhere projectdirx(args...) = joinpath($projectdir_static, args...)

for dir_type ∈ ("data", "src", "plots", "scripts", "papers")
    function_name = Symbol(dir_type * "dirx")
    @everywhere @eval begin
        $function_name(args...) = projectdirx($dir_type, args...)
    end
end

@everywhere include(srcdirx("Evolution.jl"))
@everywhere include(srcdirx("FitnessFunctions.jl"))
@everywhere include(srcdirx("DynamicalClustering.jl"))

@everywhere all_experiments = ["RepeatedEvolution_MutualInh_Dyn"]

for exp_name in all_experiments

    @everywhere include(srcdirx("ExperimentSetups/" * $exp_name * ".jl"))

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
    summaryd["Average N non-terminated"] = mean(map(et->count(x->x!=ReturnCode.Terminated ,et.retcodes) / length(et.retcodes),sim))

    @tag!(summaryd)

    safesave(datadirx("exp_summaries",exp_name * "_Summary.jld2"), summaryd)

    ########################################

    @everywhere end_networks = map(et->et.traversed_networks[end],sim);
    @everywhere end_networks_t2s = map(et->et.traversed_t2s[end],sim);

    n_networks = length(end_networks)

    end_networks_dyn_v = pmap(nt->get_rel_dyn_vector(nt[1],nt[2],n_steps,save_id),zip(end_networks,end_networks_t2s));

    end_X = reduce(hcat,end_networks_dyn_v)

    dmat_m = pairwise(d_metric,end_X,dims = 2)

    ########################################

    fundamental_networks_dyn_v = pmap(nt->get_rel_dyn_vector(nt[1],nt[2],n_steps,save_id),zip(fundamental_networks,fundamental_networks_t2s));

    fund_X = reduce(hcat,fundamental_networks_dyn_v)

    fund_dmat_m = pairwise(d_metric,end_X,fund_X,dims = 2)

    ########################################

    fulld = Dict{String, Any}()

    fulld["fitness_traj"] = map(et->et.fitness_trajectory,sim)
    fulld["t2s_traj"] = map(et->et.traversed_t2s,sim)
    fulld["geno_traj"] = map(et->reduce(hcat,map(x->vec(x),et.traversed_networks)),sim)
    fulld["retcodes"] = map(et->map(x-> x == ReturnCode.Terminated ? 1 : 0,et.retcodes),sim)
    fulld["dmat"] = dmat_m
    fulld["fund_dmat"] = fund_dmat_m
 
    @tag!(fulld)

    safesave(datadirx("exp_raw",exp_name * "_RawData.jld2"), fulld)

end
