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

cluster_calc = false

# if !cluster_calc
#     @quickactivate "GRNEvoContingency"
# end

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

@everywhere include(srcdirx("GRNEvoContingency.jl"))

all_experiments = ["DeNovoStripe/RE_Minimal_Inhibiting_DeNovo_test"]

for exp_name in all_experiments

    @everywhere include(srcdirx("ExperimentSetups/" * $exp_name * ".jl"))

    sim = pmap(worker-> SSWM_Evolution(start_network,grn_parameters,β,max_gen,tolerance,fitness_function,mutate_function),[worker for worker in 1:n_trials])

    n_traj = length(findall(x->all(x.converged),sim))

    summaryd = Dict{String, Any}()

    summaryd["Total traj simulated"] = n_trials
    summaryd["Total traj converged"] = n_traj
    summaryd["ConvergenceRate"] = n_traj/n_trials
    summaryd["N unique workers"] = length(unique(map(et->et.worker_id,sim)))
    summaryd["Average N non-terminated"] = mean(map(et->count(x->x!=ReturnCode.Terminated,et.retcodes) / length(et.retcodes),sim))

    @tag!(summaryd)

    safesave(datadirx("exp_summaries",exp_name * "_Summary.jld2"), summaryd)

    ########################################

    conv = map(x->x.converged,sim)

    print(sum(conv))

    if sum(conv) != 0

        fitness_traj_conv = map(et->et.fitness_trajectory,sim[conv]);

        end_networks = map(et->et.traversed_networks[end],sim[conv]);
        end_networks_t2s = map(et->et.traversed_t2s[end],sim[conv]);

        # stripe_achieved = map(x->minimum(findall(x->x[1] == 0.,unique(x))),fitness_traj_conv);

        stripe_achieved = map(ft->minimum(findall(x->x[1] == 0.,ft)),fitness_traj_conv);

        first_stripe_networks = [et.traversed_networks[id] for (et,id) in zip(sim[conv],stripe_achieved)];

        ########################################

        sendto(workers(), end_networks=end_networks)
        sendto(workers(), end_networks_t2s=end_networks_t2s)

        n_networks = length(end_networks)

        end_networks_dyn_cell = pmap(nt->get_rel_dyn_vector(nt[1],nt[2],n_steps,save_id),zip(end_networks,end_networks_t2s));
        end_networks_dyn_av = pmap(nt->get_av_dyn_vector(nt[1],nt[2],n_steps,n_segments),zip(end_networks,end_networks_t2s));

        print(length(end_networks_dyn_cell))

        end_X_cell = reduce(hcat,end_networks_dyn_cell)
        end_X_av = reduce(hcat,end_networks_dyn_av)

        # dmat_m_cell = pairwise(d_metric,end_X_cell,dims = 2)
        # dmat_m_av = pairwise(d_metric,end_X_av,dims = 2)


        ########################################

        sendto(workers(), first_stripe_networks=first_stripe_networks)

        min_end_networks = pmap(n->find_minimal_network(vec(n)[1:10],grn_parameters,DefaultGRNSolver(),fitness_function),end_networks);
        min_first_stripe_networks = pmap(n->find_minimal_network(vec(n)[1:10],grn_parameters,DefaultGRNSolver(),fitness_function),first_stripe_networks);

        ########################################

        fulld = Dict{String, Any}()

        fulld["fitness_traj"] = map(et->et.fitness_trajectory,sim)
        fulld["wait_times"] = map(et->et.wait_times,sim)
        fulld["t2s_traj"] = map(et->et.traversed_t2s,sim)
        fulld["geno_traj"] = map(et->reduce(hcat,map(x->vec(x),et.traversed_networks)),sim)
        fulld["retcodes"] = map(et->map(x-> x == ReturnCode.Terminated ? 1 : 0,et.retcodes),sim)
        fulld["mut_choices"] = map(et->et.mut_choices,sim)
        fulld["mut_type"] = map(et->et.mut_type,sim)
        fulld["mut_sizes"] = map(et->et.mut_sizes,sim)

        fulld["dmat_X_cell"] = end_X_cell
        fulld["dmat_X_av"] = end_X_av

        fulld["min_end_networks"] = min_end_networks
        fulld["min_fs_networks"] = min_first_stripe_networks


        fulld["fs_id"] = stripe_achieved
        fulld["converged"]  = conv
    
        @tag!(fulld)

        safesave(datadirx("exp_raw",exp_name * "_RawData.jld2"), fulld)

    else
        print("No converged trajectories. Aborting further analysis.")
    end

end
