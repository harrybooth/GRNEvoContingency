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
@everywhere include(srcdirx("FitnessLandscapes.jl"))
@everywhere include(srcdirx("DynamicalClustering.jl"))
@everywhere include(srcdirx("MinimalNetworks.jl"))

@everywhere all_experiments = ["DeNovoStripe/RE_Minimal_Inhibiting_DeNovo"]

N_interp_points = 20

for exp_name in all_experiments

    @everywhere include(srcdirx("ExperimentSetups/" * $exp_name * ".jl"))

    @everywhere development = DefaultGRNSolver()

    end_networks = load(srcdirx("ExperimentSetups/DeNovoStripe/all_end_networks.jld2"))["data"]

    sendto(workers(), end_networks=end_networks)

    n_networks = length(end_networks)

    dmat_id = [(i,j) for i in 1:n_networks, j in 1:n_networks if i>j]

    dmat = pmap(id->instability_lean(end_networks[id[1]],end_networks[id[2]],N_interp_points,grn_parameters,development,fitness_function),dmat_id)

    dmat_m = zeros(n_networks,n_networks)

    for (n,id) in enumerate(dmat_id)
        dmat_m[id...] = dmat[n]
    end

    for id in [(i,j) for i in 1:n_networks, j in 1:n_networks if i<j]
        dmat_m[id...] = dmat_m[id[2],id[1]]
    end

    ########################################

    fulld = Dict{String, Any}()

    fulld["dmat_lmc"] = dmat
 
    @tag!(fulld)

    safesave(datadirx("exp_raw",exp_name * "_LMC.jld2"), fulld)

end
