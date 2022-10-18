using DrWatson

@quickactivate "GRNEvoContingency"


using DifferentialEquations
using Random
using Parameters
using StatsBase
using Printf
using Distributed

using Base.Threads
using Base.Threads: @spawn

@everywhere include(srcdir("TissueModel_ND.jl"))
@everywhere include(srcdir("NetworkTopologies.jl"))

# Define viable mutations

viable_mutations = ones(Ng,Ng+1)

viable_mutations[2,4] = 0.
viable_mutations[3,4] = 0.

# Define randomized network generators

function random_feed_forward(prob,i,repeat)
    Random.seed!(i)
    remake(prob,p = ((0.9995 .^ rand(0:10000,Ng,Ng+1)) .* 10 .* rand(Ng,Ng+1) .* w_feed_forward,prob.p[2:end]...))
end

function random_mutual_inh(prob,i,repeat)
    Random.seed!(i)
    remake(prob,p = ((0.9995 .^ rand(0:10000,Ng,Ng+1)) .* 10 .* rand(Ng,Ng+1) .* w_mutual_inh,prob.p[2:end]...))
end

function random_frozen_osc(prob,i,repeat)
    Random.seed!(i)
    remake(prob,p = ((0.9995 .^ rand(0:10000,Ng,Ng+1)) .* 10 .* rand(Ng,Ng+1) .* w_frozen_osc,prob.p[2:end]...))
end

function random_overlap_dom(prob,i,repeat)
    Random.seed!(i)
    remake(prob,p = ((0.9995 .^ rand(0:10000,Ng,Ng+1)) .* 10 .* rand(Ng,Ng+1) .* w_overlap_dom,prob.p[2:end]...))
end

function random_bistable(prob,i,repeat)
    Random.seed!(i)
    remake(prob,p = ((0.9995 .^ rand(0:10000,Ng,Ng+1)) .* 10 .* rand(Ng,Ng+1) .* w_bistable,prob.p[2:end]...))
end

function random_classical(prob,i,repeat)
    Random.seed!(i)
    remake(prob,p = ((0.9995 .^ rand(0:10000,Ng,Ng+1)) .* 10 .* rand(Ng,Ng+1) .* w_classical,prob.p[2:end]...))
end

function random_network_bs(prob,i,repeat)
    Random.seed!(i)
    remake(prob,p = ((0.9995 .^ rand(0:10000,Ng,Ng+1)) .* 10 .* rand(Ng,Ng+1) .* viable_mutations,prob.p[2:end]...))
end

function random_network(prob,i,repeat)
    Random.seed!(i)
    remake(prob,p = (10 .* rand(Ng,Ng+1) .* viable_mutations,prob.p[2:end]...))
end

problem_types = Dict("feed_forward" => random_feed_forward,
                    "mutual_inh" => random_mutual_inh,
                    "frozen_osc" => random_frozen_osc,
                    "overlap_dom" => random_overlap_dom,
                    "bistable"=> random_bistable,
                    "classical" => random_classical,
                    "random" => random_network,
                    "random_bs" => random_network_bs);

# Define simulation output

function get_outputs(prob_ode,alg,remake_func,i,ss_abstol,ss_reltol)
    new_prob = remake_func(prob_ode,i,false)
    sol = solve(new_prob,alg, isoutofdomain=(u,p,t) -> any(x -> x < 0, u), callback = TerminateSteadyState(ss_abstol,ss_reltol),maxiters = 1e6 + 1, verbose = false, save_everystep = false)
    return (sol.retcode,sol.destats.naccept,sol.destats.nreject)
end

# Test network function

function test_networks(n_traj,method = "random",algo = "Tsit5",ss_abstol=1e-5,ss_reltol=1e-3)

    w = zeros(Ng,Ng+1)

    λg = 0.05 .* ones(Ng);

    p = (w,λg)

    g0 = init_gene_regulation_1d(0.01);

    tspan = (0,Inf)

    prob_ode = ODEProblem(gene_regulation_1d!,g0,tspan,p);

    sim = fill((:Success,1,1),n_traj)

    if algo == "Tsit5-Rosenbrock23"
        alg_choice = AutoTsit5(Rosenbrock23())
    elseif algo == "Tsit5-Radau"
        alg_choice = AutoTsit5(RadauIIA5())
    elseif algo == "Tsit5"
        alg_choice = Tsit5()
    else
        throw(DomainError(algo,"Wrong solver spec"))
    end

    @sync for i in 1:n_traj
        @spawn sim[i] = get_outputs(prob_ode,alg_choice,problem_types[method],i,ss_abstol,ss_reltol)
    end

    n_2 = count(x->x[2] + x[3] < 1e2,sim)
    n_3 = count(x->x[2] + x[3] < 1e3,sim)
    n_4 = count(x->x[2] + x[3] < 1e4,sim)
    n_5 = count(x->x[2] + x[3] < 1e5,sim)
    n_6 = count(x->x[2] + x[3] < 1e6,sim)

    n_unstable = count(x->x[1]==:Unstable,sim)
    n_maxiter = count(x->x[1]==:MaxIters,sim)

    med_it = median(map(x->x[2] + x[3] ,sim))

    return sim,n_2,n_3,n_4,n_5,n_6,n_unstable,n_maxiter,med_it
end

# Make simulation : https://juliadynamics.github.io/DrWatson.jl/dev/workflow/

function makesim(d::Dict)
    
    @unpack n_traj, method,algo,ss_abstol,ss_reltol = d

    r,n_2,n_3,n_4,n_5,n_6,n_unstable,n_maxiter,med_it = test_networks(n_traj,method,algo,ss_abstol,ss_reltol)

    fulld = copy(d)

    fulld["it<1e2"] = n_2
    fulld["it<1e3"] = n_3
    fulld["it<1e4"] = n_4
    fulld["it<1e5"] = n_5
    fulld["it<1e6"] = n_6
    fulld["median_it"] = med_it

    fulld["n_maxiters"] = n_maxiter
    fulld["n_unstable"] = n_unstable

    fulld["raw_data"] = r
    fulld["raw_data_fields"] = ["sol.retcode","sol.destats.naccept","sol.destats.nreject"]

    return fulld
end

# Run

n_traj_r = 100
n_traj_o = 100

test_specification = Dict("n_traj" =>[n_traj_r,n_traj_r,n_traj_o,n_traj_o,n_traj_o,n_traj_o,n_traj_o,n_traj_o],
                          "method"=>["random","random_bs","classical","bistable","overlap_dom","frozen_osc","mutual_inh","feed_forward"],
                          "algo"=>"Tsit5-Radau",
                          "ss_abstol" => 1e-5,
                          "ss_reltol" => 1e-3)

all_tests = dict_list(test_specification);

for (i,d) in enumerate(all_tests)
    f = makesim(d)
    safesave(datadir("sims/termination_testing", savename(d, "jld2")), f)
end