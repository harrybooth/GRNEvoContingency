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

# http://jpfairbanks.com/2017/12/27/running-julia-on-slurm-cluster/
# 
# This code is expected to be run from an sbatch script after a module load julia command has been run.
# It starts the remote processes with srun within an allocation specified in the sbatch script.

# https://github.com/JuliaLang/julia/pull/30174
# https://github.com/JuliaLang/julia/issues/31953

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

@everywhere all_experiments = ["FindNetworks_HalfStripeLeft"]

for exp_name in all_experiments

    @everywhere include(srcdirx("ExperimentSetups/" * $exp_name * ".jl"))

    fulld = Dict{String, Any}()

    for network_topology in networks_to_search

        summaryd = Dict{String, Any}()

        @everywhere begin

            start_top = network_topology_dict[$network_topology]
            
            viable_mutations = Int64.(start_top .!= 0)

            mutation_weights = findall(viable_mutations .> 0)

            n_sample_func() = rand(Binomial(length(mutation_weights),mut_prob))

            mutation_op = MutationOperator(Normal,(μ = 0.0,σ = noise_cv),n_sample_func,deletion_prob,max_w,mutation_weights)

            mutate_function = i -> noise_no_additions(i,mutation_op);
        end

        evo_trace = SSWM_Evolution((0.9995 .^ rand(0:10000,Ng,Ng+1)) .* max_w .* rand(Ng,Ng+1) .* start_top,grn_parameters,β,max_gen,tolerance,fitness_function,mutate_function)

        evo_trace.converged = false

        sim = fill(evo_trace,n_networks_required)

        sim_id = findall(x->!x.converged,sim)

        n_trials =  0

        while length(sim_id) != 0

            n_trials += length(sim_id)

            sim = pmap(sim_i -> sim_i.converged ? sim_i : SSWM_Evolution((0.9995 .^ rand(0:10000,Ng,Ng+1)) .* max_w .* rand(Ng,Ng+1) .* start_top,grn_parameters,β,max_gen,tolerance,fitness_function,mutate_function),sim)

            sim_id = findall(x->!x.converged,sim)

        end

        fulld[network_topology * "_networks"] = map(et->et.traversed_networks[end],sim);
        fulld[network_topology * "_fitness"] = map(et->et.fitness_trajectory[end],sim);
        fulld[network_topology * "_t2s"] = map(et->et.traversed_t2s[end],sim);

        summaryd["Total traj simulated"] = n_trials
        summaryd["Total traj converged"] = n_networks_required
        summaryd["ConvergenceRate"] = n_networks_required/n_trials
        summaryd["N unique workers"] = length(unique(map(et->et.worker_id,sim)))
        summaryd["Average N non-terminated"] = mean(map(et->count(x->x!=ReturnCode.Terminated,et.retcodes) / length(et.retcodes),sim))

        @tag!(summaryd)

        safesave(datadirx("exp_summaries",exp_name * "_" * network_topology * "_Summary.jld2"), summaryd)
    end

    @tag!(fulld)

    safesave(datadirx("exp_raw",exp_name * "_RawData.jld2"), fulld)

end
