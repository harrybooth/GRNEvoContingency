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

@everywhere all_experiments = ["RepeatedEvolution_Bistable_Dyn"]

function last_networks(run_id,fitness_route,geno_route) 

    unique_fitness_routes = unique(fitness_route)

    result = []

    for fb in unique_fitness_routes

        fb_ids = findall(x->x == fb,fitness_route)

        last_id = maximum(fb_ids)

        push!(result,(run_id,fb,geno_route[:,last_id]))
    end

    return result

end

for exp_name in all_experiments

    @everywhere include(srcdirx("ExperimentSetups/" * $exp_name * ".jl"))

    exp_data = load(datadirx("exp_raw",exp_name * "_RawData.jld2"))

    ############

    fitness_traj = copy(exp_data["fitness_traj"])
    geno_traj = copy(exp_data["geno_traj"])
    
    all_unique_fitness = unique(reduce(vcat,map(x->unique(x),fitness_traj)))
    
    initial_fitness = fitness_traj[1][1]
    
    n_bin = 10
    
    hist_edges = zeros(n_bin+1)
    
    hist_edges[1] = initial_fitness
    
    hist_edges[2:n_bin] .= LinRange(initial_fitness+eps(),0.9,n_bin-1) |> collect
    
    hist_edges[n_bin+1] = 1.
    
    n_fit_bin = length(hist_edges) - 1
    
    h_fitness = fit(Histogram, all_unique_fitness, hist_edges; closed = :left) 
    
    fitness_routes = map(traj->map(f->StatsBase.binindex(h_fitness, f),unique(traj)),fitness_traj);
    
    ln_fr = map(x -> last_networks(x[1],x[2][1],x[2][2]), enumerate(zip(fitness_routes,geno_traj)));
    
    ln_frv = reduce(vcat,ln_fr)

    ###############

    result =  Dict()

    for fb_n in 2:n_bin-1

        fb_networks = filter(x->x[2] == fb_n, ln_frv)

        eval_networks = map(x->x[3],fb_networks)

        sendto(workers(), eval_networks=eval_networks)

        run_ids = map(x->x[1],fb_networks)

        n_sample = length(eval_networks)

        dmat_id = [(i,j) for i in 1:n_sample, j in 1:n_sample if i>j]

        dmat = pmap(id-> instability_lean(eval_networks[id[1]],eval_networks[id[2]],N_interp_points,grn_parameters,development,fitness_function),dmat_id)

        dmat_m = zeros(n_sample,n_sample)

        for (n,id) in enumerate(dmat_id)
            dmat_m[id...] = dmat[n]
        end

        for id in [(i,j) for i in 1:n_networks, j in 1:n_networks if i<j]
            dmat_m[id...] = dmat_m[id[2],id[1]]
        end

        result[fb_n] = (dmat_m,run_ids,eval_networks)

        print("Status: completed " * string(fb_n) * " \n")

    end

    ########################################

    fulld = Dict{String, Any}()

    fulld["dmat_lmc_binned"] = result
 
    @tag!(fulld)

    safesave(datadirx("exp_raw",exp_name * "_LMC.jld2"), fulld)

end
