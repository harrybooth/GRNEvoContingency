using Pkg

Pkg.activate("../..")

projectdir_static = dirname(Base.active_project())

projectdirx(args...) = joinpath(projectdir_static, args...)

# Generate functions to access the path of default subdirectories.
for dir_type ∈ ("data", "src", "plots", "scripts", "papers")
    function_name = Symbol(dir_type * "dirx")
    @eval begin
        $function_name(args...) = projectdirx($dir_type, args...)
    end
end

using DrWatson

@quickactivate "GRNEvoContingencyAnalysis"

projectname()

include(srcdirx("GRNEvoContingency.jl"))
include(srcdirx("Analysis.jl"))

using JLD2
using CairoMakie
using ColorSchemes
using BenchmarkTools
using Distances
using StatsPlots
using Clustering
using MultivariateStats
using HypothesisTests
using PyCall
using DecisionTree
using DataFrames
using Combinatorics
using GraphMakie
using Graphs
using NetworkLayout
using LinearAlgebra
using UMAP
using PyCallJLD2

xgboost = pyimport("xgboost");

all_conv = []
all_geno_traj = []
all_fitness_traj = []
all_min_end_networks = []
all_min_fs_networks = []
all_wait_times = []
all_mut_types = []
all_mut_sizes = []
all_mut_choices = []
all_mat_av_X = []

exp_versions = ["N1","N2","N3","N4","N5","N6"]

for ver in exp_versions

    data = load(datadirx("exp_pro/DeNovoStripe/" * exp_name * "_RawData_" * ver * ".jld2"));

    conv = data["converged"]

    push!(all_conv,conv)
    push!(all_geno_traj,data["geno_traj"][conv])
    push!(all_fitness_traj,data["fitness_traj"][conv])
    push!(all_min_end_networks,data["min_end_networks"])
    push!(all_min_fs_networks,data["min_fs_networks"])
    push!(all_wait_times,data["wait_times"][conv])
    push!(all_mut_types,data["mut_type"][conv])
    push!(all_mut_sizes,data["mut_sizes"][conv])
    push!(all_mut_choices,data["mut_choices"][conv])
    push!(all_mat_av_X,data["dmat_X_av"])
end

all_conv = reduce(vcat,all_conv)
all_geno_traj = reduce(vcat,all_geno_traj)
all_fitness_traj = reduce(vcat,all_fitness_traj)
all_min_end_networks = reduce(vcat,all_min_end_networks)
all_min_fs_networks = reduce(vcat,all_min_fs_networks)
all_wait_times = reduce(vcat,all_wait_times)
all_mut_types = reduce(vcat,all_mut_types)
all_mut_sizes = reduce(vcat,all_mut_sizes)
all_mut_choices = reduce(vcat,all_mut_choices)
all_mat_av_X = reduce(hcat,all_mat_av_X);

fs_mss =  map(list_mss->select_minimal_topologies(list_mss),all_min_fs_networks)
ls_mss =  map(list_mss->select_minimal_topologies(list_mss),all_min_end_networks);

nconv = sum(all_conv)

###################################

trajectories = map(n->Trajectory(n,all_geno_traj[n],all_fitness_traj[n],all_wait_times[n] .+ 1,all_mut_choices[n],all_mut_types[n],all_mut_sizes[n],weight_names),1:nconv);

for (n,tr) in enumerate(trajectories)
    assign_minimal_subgraphs!(tr,fs_mss[n],ls_mss[n])
end

inc_metagraph, vertex_top_map,top_vertex_map,vertex_complexity_map,inclusion_matrix = create_inclusion_metagraph(trajectories);

minimal_motif_id = findall(indegree(inc_metagraph) .== 0);

minimal_motifs = reduce(hcat,[vertex_top_map[vertex_id] for vertex_id in minimal_motif_id])

for tr in trajectories
    assign_inc_vertex_ids!(tr,top_vertex_map)
    assign_inc_parents!(tr,inclusion_matrix,vertex_complexity_map,minimal_motif_id)
end

##################################

train_ttl = false
train_gtl = false

top_n = 4

predict_id = sorted_uep[1:top_n]

label_names = 1:top_n |> collect

trajectories_p = copy(trajectories)

vertex_to_predict_label = Dict(vertex=>n for (n,vertex) in enumerate(predict_id))
predict_label_to_vertex = Dict(n=>vertex for (n,vertex) in enumerate(predict_id))
predict_label_to_vertex[top_n+1] = -1
labels = map(tr->tr.inc_metagraph_vertices[tr.H0] ∈ predict_id ? vertex_to_predict_label[tr.inc_metagraph_vertices[tr.H0]] : top_n + 1,trajectories_p);

###################################

const c_types = ["c" for _ in 1:10];

null_H0_dist = [count(x->x==i,labels) for i in 1:top_n+1] ./ length(trajectories_p);

ac = "DN" * string(top_n)

if train_gtl || train_ttl
    train_id,test_id =  create_train_test_id_split(labels,0.8);
    save("train_test_ids_Other_" * ac * ".jld2", Dict("train" => train_id, "test" => test_id))
else
    train_id = load(local_nb_data * "train_test_ids_Other_" * ac * ".jld2", "train")
    test_id = load(local_nb_data * "train_test_ids_Other_" * ac * ".jld2", "test")
end;

##################################

X_train_ttl_v = [reduce(hcat,unique([vcat(features,label) for features in tr.topologies[1:tr.H0-1]])) for (label,tr) in zip(labels[train_id],trajectories_p[train_id])]
X_train_ttl = reduce(hcat,X_train_ttl_v) |> transpose |> collect

y_train_ttl = copy(Int.(X_train_ttl[:,13]))

X_train_gtl_v = [reduce(hcat,[vcat(features,label) for features in tr.geno_traj[1:tr.H0-1]]) for (label,tr) in zip(labels[train_id],trajectories_p[train_id])]
X_train_gtl  = reduce(hcat,X_train_gtl_v) |> transpose |> collect;

y_train_gtl = copy(Int.(X_train_gtl[:,13]))

#################################

X_test_ttl_v = [reduce(hcat,unique([vcat(features,label) for features in tr.topologies[1:tr.H0-1]])) for (label,tr) in zip(labels[test_id],trajectories_p[test_id])]
X_test_ttl = reduce(hcat,X_test_ttl_v) |> transpose |> collect

y_test_ttl = copy(Int.(X_test_ttl[:,13]))

X_test_gtl_v = [reduce(hcat,[vcat(features,label) for features in tr.geno_traj[1:tr.H0-1]]) for (label,tr) in zip(labels[test_id],trajectories_p[test_id])]
X_test_gtl  = reduce(hcat,X_test_gtl_v) |> transpose |> collect;

y_test_gtl = copy(Int.(X_test_gtl[:,13]));

#################################

if train_ttl

    params = Dict(
        "eta"=> 0.01,
        "objective"=>"multi:softprob",
        "num_class"=>top_n+1,
        "subsample"=> 0.5,
        "eval_metric"=>"auc"
    )

    max_boosting_rounds = 5000
    early_stop = 200

    d_train_ttl = xgboost.DMatrix(X_train_ttl[:,1:10], label=y_train_ttl .- 1,feature_names = weight_names)
    d_test_ttl = xgboost.DMatrix(X_test_ttl[:,1:10], label=y_test_ttl .- 1,feature_names = weight_names)

    model_ttl = xgboost.train(params, d_train_ttl, max_boosting_rounds, evals = [(d_test_ttl, "test")], verbose_eval=false, early_stopping_rounds=early_stop)

    ttl_prob_train = model_ttl.predict(d_train_ttl)
    ttl_prob_test = model_ttl.predict(d_test_ttl)

    y_pred_train_ttl = mapslices(x->argmax(x),ttl_prob_train,dims = 2)
    y_pred_test_ttl = mapslices(x->argmax(x),ttl_prob_test,dims = 2);

    print("Accuracy (test): " * string(sum(y_pred_test_ttl .== y_test_ttl) / length(y_test_ttl)))
    print("\n")

    all_bs_test = []

    for (n,p_row) in enumerate(eachrow(ttl_prob_test))
        push!(all_bs_test,brier_score(p_row,y_test_ttl[n]))
    end

    print("Brier score (test): " * string(mean(all_bs_test)))
    print("\n")

    all_bs_train = []

    for (n,p_row) in  enumerate(eachrow(ttl_prob_train))
        push!(all_bs_train,brier_score(p_row,y_train_ttl[n]))
    end

    print("Brier score (train): " * string(mean(all_bs_train)))
    print("\n")
    
    model_ttl.save_model("ModelTTL_Other_" * ac * ".json")

else
    model_ttl = xgboost.Booster()

    d_train_ttl = xgboost.DMatrix(X_train_ttl[:,1:10], label=y_train_ttl .- 1,feature_names = weight_names)
    d_test_ttl = xgboost.DMatrix(X_test_ttl[:,1:10], label=y_test_ttl .- 1,feature_names = weight_names)

    model_ttl.load_model(local_nb_data * "ModelTTL_Other_" * ac * ".json")

    ttl_prob_train = model_ttl.predict(d_train_ttl)
    ttl_prob_test = model_ttl.predict(d_test_ttl)

    y_pred_train_ttl = mapslices(x->argmax(x),ttl_prob_train,dims = 2)
    y_pred_test_ttl = mapslices(x->argmax(x),ttl_prob_test,dims = 2);
end;


if train_gtl

    params = Dict(
        "eta"=> 0.01,
        "objective"=>"multi:softprob",
        "num_class"=>top_n+1,
        "subsample"=> 0.5,
        "eval_metric"=>"auc"
    )

    max_boosting_rounds = 5000
    early_stop = 200

    d_train_gtl = xgboost.DMatrix(X_train_gtl[:,1:10], label=y_train_gtl .- 1,feature_names = weight_names)
    d_test_gtl = xgboost.DMatrix(X_test_gtl[:,1:10], label=y_test_gtl .- 1,feature_names = weight_names)

    model_gtl = xgboost.train(params, d_train_gtl, max_boosting_rounds, evals = [(d_test_gtl, "test")], verbose_eval=false, early_stopping_rounds=early_stop)

    gtl_prob_train = model_gtl.predict(d_train_gtl)
    gtl_prob_test = model_gtl.predict(d_test_gtl)

    y_pred_train_gtl = mapslices(x->argmax(x),gtl_prob_train,dims = 2)
    y_pred_test_gtl = mapslices(x->argmax(x),gtl_prob_test,dims = 2);

    print("Accuracy (test): " * string(sum(y_pred_test_gtl .== y_test_gtl) / length(y_test_gtl)))
    print("\n")

    all_bs_test = []

    for (n,p_row) in  enumerate(eachrow(gtl_prob_test))
        push!(all_bs_test,brier_score(p_row,y_test_gtl[n]))
    end

    print("Brier score (test): " * string(mean(all_bs_test)))
    print("\n")

    all_bs_train = []

    for (n,p_row) in  enumerate(eachrow(gtl_prob_train))
        push!(all_bs_train,brier_score(p_row,y_train_gtl[n]))
    end

    print("Brier score (train): " * string(mean(all_bs_train)))
    print("\n")

    model_gtl.save_model("ModelGTL_Other_" * ac * ".json")

else
    model_gtl = xgboost.Booster()

    d_train_gtl = xgboost.DMatrix(X_train_gtl[:,1:10], label=y_train_gtl .- 1,feature_names = weight_names)
    d_test_gtl = xgboost.DMatrix(X_test_gtl[:,1:10], label=y_test_gtl .- 1,feature_names = weight_names)

    model_gtl.load_model(local_nb_data * "ModelGTL_Other_" * ac * ".json")

    gtl_prob_train = model_gtl.predict(d_train_gtl)
    gtl_prob_test = model_gtl.predict(d_test_gtl)

    y_pred_train_gtl = mapslices(x->argmax(x),gtl_prob_train,dims = 2)
    y_pred_test_gtl = mapslices(x->argmax(x),gtl_prob_test,dims = 2);

end;

#########################

include(srcdirx("TrajectoryAnalysis.jl"))

for tr in trajectories_p
    assign_predictions!(tr,model_ttl,:tt,predict_label_to_vertex)
    assign_predictions!(tr,model_gtl,:gt,predict_label_to_vertex)
end

for (tr,label) in zip(trajectories_p,labels)
    assign_tt_other_prediction_errors!(tr,predict_label_to_vertex[label],predict_id)
    assign_gt_other_prediction_errors!(tr,predict_label_to_vertex[label],predict_id)
end

for tr in trajectories_p
    assign_weight_edits!(tr)
end

