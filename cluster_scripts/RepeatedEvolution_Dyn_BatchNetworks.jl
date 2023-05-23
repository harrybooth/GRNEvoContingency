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
end

@everywhere begin
    using DrWatson
    using JLD2
    using Printf
    using Base.Threads
    using Base.Threads: @spawn
    using ParallelDataTransfer
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

@everywhere all_experiments = ["RepeatedEvolution_FrozenOsc_Dyn_Batch"]

for exp_name in all_experiments

    @everywhere include(srcdirx("ExperimentSetups/RepeatedEvolution/" * $exp_name * ".jl"))

    shift = 4

    for choice in 1+shift:n_test_networks+shift

        start_network = start_networks_dict[topology_choice * "_networks"][choice]

        sendto(workers(), start_network=start_network)

        sim = []

        n_trials =  0

        while length(sim) != n_traj

            sim_temp = pmap(worker-> SSWM_Evolution(start_network,grn_parameters,β,max_gen,tolerance,fitness_function,mutate_function),[worker for worker in 1:n_tasks])

            n_trials += n_tasks

            sim_temp_id = findall(x->x.converged,sim_temp)

            for id in sim_temp_id
                if length(sim) < n_traj
                    push!(sim,sim_temp[id])
                end
            end
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

        end_networks = map(et->et.traversed_networks[end],sim);
        end_networks_t2s = map(et->et.traversed_t2s[end],sim);

        sendto(workers(), end_networks=end_networks)
        sendto(workers(), end_networks_t2s=end_networks_t2s)

        n_networks = length(end_networks)

        end_networks_dyn_cell = pmap(nt->get_rel_dyn_vector(nt[1],nt[2],n_steps,save_id),zip(end_networks,end_networks_t2s));
        end_networks_dyn_av = pmap(nt->get_av_dyn_vector(nt[1],nt[2],n_steps,n_segments),zip(end_networks,end_networks_t2s));

        end_X_cell = reduce(hcat,end_networks_dyn_cell)
        end_X_av = reduce(hcat,end_networks_dyn_av)

        dmat_m_cell = pairwise(d_metric,end_X_cell,dims = 2)
        dmat_m_av = pairwise(d_metric,end_X_av,dims = 2)

        ########################################

        fundamental_networks_dyn_cell = pmap(nt->get_rel_dyn_vector(nt[1],nt[2],n_steps,save_id),zip(fundamental_networks,fundamental_networks_t2s));
        fundamental_networks_dyn_av = pmap(nt->get_av_dyn_vector(nt[1],nt[2],n_steps,n_segments),zip(fundamental_networks,fundamental_networks_t2s));

        fund_X_cell = reduce(hcat,fundamental_networks_dyn_cell)
        fund_X_av = reduce(hcat,fundamental_networks_dyn_av)

        fund_m_cell = pairwise(d_metric,fund_X_cell,dims = 2)
        fund_m_av = pairwise(d_metric,fund_X_av,dims = 2)

        fund_dmat_m_cell = pairwise(d_metric,end_X_cell,fund_X_cell,dims = 2)
        fund_dmat_m_av = pairwise(d_metric,end_X_av,fund_X_av,dims = 2)

        ########################################

        fulld = Dict{String, Any}()

        fulld["fitness_traj"] = map(et->et.fitness_trajectory,sim)
        fulld["t2s_traj"] = map(et->et.traversed_t2s,sim)
        fulld["geno_traj"] = map(et->reduce(hcat,map(x->vec(x),et.traversed_networks)),sim)
        fulld["retcodes"] = map(et->map(x-> x == ReturnCode.Terminated ? 1 : 0,et.retcodes),sim)

        fulld["dmat_cell"] = dmat_m_cell
        fulld["dmat_X_cell"] = end_X_cell
        fulld["fund_cell"] = fund_m_cell
        fulld["fund_X_cell"] = fund_X_cell
        fulld["fund_dmat_cell"] = fund_dmat_m_cell

        fulld["dmat_av"] = dmat_m_av
        fulld["dmat_X_av"] = end_X_av
        fulld["fund_av"] = fund_m_av
        fulld["fund_X_av"] = fund_X_av
        fulld["fund_dmat_av"] = fund_dmat_m_av
    
        @tag!(fulld)

        safesave(datadirx("exp_raw",exp_name * "_RawData_" * string(choice) * ".jld2"), fulld)
    end

end
