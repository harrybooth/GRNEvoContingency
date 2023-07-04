# http://jpfairbanks.com/2017/12/27/running-julia-on-slurm-cluster/
# 
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

@everywhere all_experiments = ["Experiment_1/RE_Minimal_Activating_Single","Experiment_1/RE_Minimal_Inhibiting_Single"]

for exp_name in all_experiments

    @everywhere include(srcdirx("ExperimentSetups/" * $exp_name * ".jl"))

    fulld = Dict{String, Any}()

    network_list = []

    conv_rate_list = []

    while conv_rate <= conv_rate_required

        random_network = (0.9995 .^ rand(0:param_N,Ng,Ng+1)) .* max_w .* rand(Ng,Ng+1) .* start_top

        sim = pmap(worker-> SSWM_Evolution(random_network,grn_parameters,β,max_gen_search,tolerance,fitness_function,mutate_function),[worker for worker in 1:n_trials_search])

        n_traj = length(findall(x->all(x.converged),sim))
    
        conv_rate = n_traj/n_trials_search

        push!(network_list,random_network)
        push!(conv_rate_list,conv_rate)
    end

    fulld["networks"] = network_list
    fulld["conv_rate"] = conv_rate_list

    @tag!(fulld)

    safesave(datadirx("exp_raw",exp_name * "_FindNetworks.jld2"), fulld)

end
