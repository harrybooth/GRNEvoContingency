using Distributed
using ClusterManagers
using DrWatson

@quickactivate GRNEvoContingency

const data_dir = datadir()
const scr_dir = scrdir()
const proj_dir = projectdir()

# http://jpfairbanks.com/2017/12/27/running-julia-on-slurm-cluster/
# 
# This code is expected to be run from an sbatch script after a module load julia command has been run.
# It starts the remote processes with srun within an allocation specified in the sbatch script.

const n_cores = 100

addprocs(SlurmManager(n_cores))

@everywhere using Pkg

@everywhere Pkg.activate("..")

@everywhere Pkg.instantiate()

@everywhere Pkg.precompile()

@everywhere begin
    using DrWatson
    using JLD2
    using Printf
    using Base.Threads
    using Base.Threads: @spawn
end

@everywhere include(srcdir("Evolution.jl"))
@everywhere include(srcdir("FitnessFunctions.jl"))

@everywhere all_experiments = ["RepeatedEvolution_Experiment_1"]

for exp_name in all_experiments

    @everywhere include(srcdir("ExperimentSetups/" * $exp_name * ".jl"))

    evo_trace = SSWM_Evolution(start_network,grn_parameters,β,max_gen,tolerance,fitness_function,mutate_function)

    evo_trace.converged = false

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

    summaryd["Total traj simulated"] = n_trials
    summaryd["Total traj converged"] = n_traj
    summaryd["ConvergenceRate"] = n_traj/n_trials
    summaryd["N unique workers"] = length(unique(map(et->et.worker_id,sim)))
    summaryd["Average N non-terminated"] = mean(map(et->count(x->x!=:Terminated,et.retcodes) / length(et.retcodes),sim))

    @tag!(summaryd)

    safesave(datadir("exp_summaries",exp_name * "_Summary.jld2"), summaryd)

end
