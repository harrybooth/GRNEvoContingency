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

@everywhere all_experiments = ["MSelection/MSelection_MutualInh_Grad_Left_Stripe_1","MSelection/MSelection_MutualInh_Grad_Left_Stripe_2","MSelection/MSelection_MutualInh_Grad_Left_Stripe_3","MSelection/MSelection_MutualInh_Grad_Left_Stripe_4","MSelection/MSelection_MutualInh_Grad_Left_Stripe_5","MSelection/MSelection_MutualInh_Grad_Left_Stripe_6"]

# @everywhere all_experiments = ["RightHanded/RepeatedEvolution_Bistable_Mut_1_RH","RightHanded/RepeatedEvolution_Bistable_Mut_2_RH","RightHanded/RepeatedEvolution_Bistable_Mut_3_RH","RightHanded/RepeatedEvolution_Bistable_Mut_4_RH","RightHanded/RepeatedEvolution_Bistable_Mut_5_RH","RightHanded/RepeatedEvolution_Bistable_Mut_6_RH","RightHanded/RepeatedEvolution_Bistable_Mut_7_RH","RightHanded/RepeatedEvolution_Bistable_Mut_8_RH"]
# @everywhere all_experiments = ["RightHanded/RepeatedEvolution_FeedForward_Mut_1_RH","RightHanded/RepeatedEvolution_FeedForward_Mut_2_RH","RightHanded/RepeatedEvolution_FeedForward_Mut_3_RH","RightHanded/RepeatedEvolution_FeedForward_Mut_4_RH","RightHanded/RepeatedEvolution_FeedForward_Mut_5_RH","RightHanded/RepeatedEvolution_FeedForward_Mut_6_RH","RightHanded/RepeatedEvolution_FeedForward_Mut_7_RH","RightHanded/RepeatedEvolution_FeedForward_Mut_8_RH"]

for exp_name in all_experiments

    @everywhere include(srcdirx("ExperimentSetups/MutationTesting/" * $exp_name * ".jl"))

    sim = pmap(worker->SSWM_MSelection(start_network,grn_parameters,β,max_gen,tolerance_list,fitness_function_list,mutate_function),[worker for worker in 1:n_trials])

    n_traj = length(findall(x->all(x.converged),sim))

    summaryd = Dict{String, Any}()

    summaryd["Total traj simulated"] = n_trials
    summaryd["Total traj converged"] = n_traj
    summaryd["ConvergenceRate"] = n_traj/n_trials
    summaryd["N unique workers"] = length(unique(map(et->et.worker_id,sim)))
    summaryd["Average N non-terminated"] = mean(map(et->count(x->x!=ReturnCode.Terminated ,et.retcodes) / length(et.retcodes),sim))

    @tag!(summaryd)

    safesave(datadirx("exp_summaries",exp_name * "_Summary.jld2"), summaryd)

    ########################################

    conv_time = map(x->length(x.fitness_trajectory),sim)

    cum_conv = [sum(conv_time .< i)/n_trials for i in 1:max_gen]

    ########################################

    fulld = Dict{String, Any}()

    fulld["convergence_times"] = conv_time
    fulld["cumulative_convergence_time"] = cum_conv
    fulld["fitness_trajectories"] = map(x->x.fitness_trajectory,sim)
 
    @tag!(fulld)

    safesave(datadirx("exp_raw",exp_name * "_RawData.jld2"), fulld)

end
