#=
ml Julia/1.8.2-linux-x86_64
=#

using DrWatson
using Distributed
using ClusterManagers

@quickactivate "GRNEvoContingency"

const n_cores = 64

const n_workers = 2

addprocs(SlurmManager(n_cores), N=n_workers, topology=:master_worker, mem_per_cpu= "4G", exeflags="--project=.",time = "00:30:00",partition = "cpu")

# @everywhere using Pkg

# @everywhere pkg"activate ."

# @everywhere pkg"instantiate"

# @everywhere pkg"precompile"

@everywhere begin
    using JLD2
    using Printf
    using Base.Threads
    using Base.Threads: @spawn
end

@everywhere include(srcdir("Evolution.jl"))
@everywhere include(srcdir("FitnessFunctions.jl"))

@everywhere all_experiments = ["Experiment_1"]

for exp_name in all_experiments

    @everywhere include(srcdir("ExperimentSetups/" * $exp_name * ".jl"))

    evo_trace = SSWM_Evolution(start_network,grn_parameters,β,max_gen,tolerance,fitness_function,mutate_function)

    sim = fill(evo_trace,n_traj)

    sim_id = findall(x->!x.converged,sim)

    n_trials =  0

    while length(sim_id) != 0

        n_trials += length(sim_id)

        @sync for i in sim_id
            @spawn sim[i] = SSWM_Evolution(start_network,grn_parameters,β,max_gen,tolerance,fitness_function,mutate_function)
        end

        sim_id = findall(x->!x.converged,sim)

    end

    fulld = Dict{String, Any}()

    fulld["fitness_traj"] = map(et->et.fitness_trajectory,sim)
    fulld["geno_traj"] = map(et->reduce(hcat,map(x->vec(x),et.traversed_networks)),sim);
    fulld["retcodes"] = map(et->map(x-> x == :Terminated ? 1 : 0,et.retcodes),sim)

    @tag!(fulld)

    safesave(datadir("exp_raw",exp_name * "_RawData.jld2"), fulld)

    summaryd = Dict{String, Any}()

    summaryd["ConvergenceRate"] = n_traj/n_trials
    summaryd["N unique workers"] = length(unique(map(et->et.worker_id,sim)))
    summaryd["Average N non-terminated"] = mean(map(et->count(x->x!=:Terminated,et.retcodes) / length(et.retcodes),sim))

    @tag!(summaryd)

    safesave(datadir("exp_summaries",exp_name * "_Summary.jld2"), summaryd)

end