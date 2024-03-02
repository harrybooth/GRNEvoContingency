const mut_choice_map = reshape(1:12 |> collect,(3,4))

mutable struct Trajectory

    sim_id :: Int64

    geno_traj :: Vector{Vector{Float64}}
    topologies :: Vector{Vector{Int64}}

    n_accepted_mutants :: Int64
    acceptance_ratio :: Float64

    mutation_number :: Vector{Int64}
    stripe_indicator :: Vector{Bool}
    H0 :: Int64
    wait_times :: Vector{Int64}

    fitness_traj :: Vector{Float64}
    fitness_traj_tuple :: Any

    top_edits :: Vector{Int64}
    weight_edits :: Any
    masked_hamming_distance_H0 :: Any
    masked_hamming_distance :: Any

    initial_fitness :: Float64
    final_fitness :: Float64
    mutant_info :: Any

    inc_metagraph_vertices :: Any
    inc_metagraph_parents :: Any 
    minimal_stripe_subgraphs :: Any

    parent_inclusion_indicator :: Any

    tt_label_probabilities :: Any
    tt_label_predictions :: Any
    tt_label_entropies :: Any
    tt_prediction_error :: Any

    gt_label_probabilities :: Any
    gt_label_predictions :: Any
    gt_label_entropies :: Any
    gt_prediction_error :: Any

    mss_probabilities :: Any
    mss_predictions :: Any
    mss_entropies :: Any
    mss_prediction_error :: Any

    train_test_indicator :: Any
    train_test_indicator_mss :: Any
    embeddings :: Any

    epistasis :: Any

    tt_kl_div :: Any
    gt_kl_div :: Any
    mss_kl_div :: Any

    gt_shap :: Any
    tt_shap :: Any
    mss_shap :: Any

    other :: Any
    
end

function Trajectory(sim_id::Int64,geno_traj_m::Matrix{Float64},fitness_traj_tuple::Vector{Tuple{Float64, Float64}},wait_times,mut_types::Any,weight_names)

    geno_traj = [collect(i) for i in eachcol(geno_traj_m)]

    ######## process trajectory data

    topologies = map(w->Int.(sign.(w)),geno_traj)

    wait_times_v = vcat([1.],wait_times)

    n_accepted_mutants = length(fitness_traj_tuple)-1
    n_generated_mutants = length(fitness_traj_tuple)-1
    acceptance_ratio = n_accepted_mutants / n_generated_mutants

    mutation_number = [i for i in 0:n_accepted_mutants]
    stripe_indicator =  map(ft->ft[1] == 0,fitness_traj_tuple)
    H0 = minimum(findall(stripe_indicator))

    fitness_traj_add = map(ft->add_fitness(ft),fitness_traj_tuple)
    top_edits = compute_cumulative_edits(reduce(hcat,topologies))
    weight_edits = nothing 
    masked_hamming_distance_H0 = nothing
    masked_hamming_distance = nothing

    initial_fitness = fitness_traj_add[1]
    final_fitness = fitness_traj_add[end]

    ######## MST

    minimal_stripe_subgraphs = nothing 
    metagraph_vertices = nothing
    metagraph_parents = nothing

    ######## generate mutation data

    weight_id = get_mutation_id(geno_traj_m)
    weight_noise = get_mutant_dist(geno_traj_m)

    fitness_delta = get_fitness_delta(fitness_traj_add)

    n_weight_changes = map(m->length(m),weight_id)

    weight_id_label = map(v->join(map(x->weight_names[x],v),"|"),weight_id)

    start_fitness = fitness_traj_add[1:end-1]
    start_fitness_tuple = fitness_traj_tuple[1:end-1]
    mutant_fitness = fitness_traj_add[2:end]

    start_network = [collect(v) for v in eachcol(geno_traj_m[:,1:end-1])]

    mutant_info = [(weight_id = weight_id[n],weight_id_label = weight_id_label[n],mut_type = mut_types[n], mut_size = weight_noise[n],start_fitness = start_fitness[n],start_fitness_tuple = start_fitness_tuple[n],mutant_fitness = mutant_fitness[n],fitness_delta = fitness_delta[n],start_network = start_network[n]) for n in 1:length(weight_id)]

    ####### predictions

    parent_inclusion_indicator = nothing

    tt_label_probabilities = nothing
    tt_label_predictions = nothing
    tt_label_entropies = nothing
    tt_prediction_error = nothing

    gt_label_probabilities = nothing
    gt_label_predictions = nothing
    gt_label_entropies = nothing
    gt_prediction_error = nothing

    mss_probabilities = nothing
    mss_predictions = nothing
    mss_entropies = nothing
    mss_prediction_error = nothing

    train_test_indicator = nothing
    train_test_indicator_mss = nothing
    embeddings = nothing
    other = nothing

    ####### instantiate 

    Trajectory(sim_id,geno_traj,topologies,n_accepted_mutants,acceptance_ratio,mutation_number,stripe_indicator,H0,wait_times_v,fitness_traj_add,fitness_traj_tuple,top_edits,weight_edits,masked_hamming_distance_H0,masked_hamming_distance,initial_fitness,final_fitness,mutant_info,metagraph_vertices,metagraph_parents,minimal_stripe_subgraphs,parent_inclusion_indicator,
                                                                                                            tt_label_probabilities,tt_label_predictions,tt_label_entropies,tt_prediction_error,gt_label_probabilities,gt_label_predictions,gt_label_entropies,gt_prediction_error,
                                                                                                            mss_probabilities,mss_predictions,mss_entropies,mss_prediction_error,train_test_indicator,train_test_indicator_mss,embeddings,other)
end

function Trajectory(sim_id::Int64,geno_traj_m::Matrix{Float64},fitness_traj_tuple::Vector{Tuple{Float64, Float64}},wait_times,mut_choices::Any,mut_types::Any,mut_sizes::Any,weight_names)

    geno_traj = [collect(i) for i in eachcol(geno_traj_m)]

    ######## process trajectory data

    topologies = map(w->Int.(sign.(w)),geno_traj)

    wait_times_v = vcat([1.],wait_times)

    n_accepted_mutants = length(fitness_traj_tuple)-1
    n_generated_mutants = length(fitness_traj_tuple)-1
    acceptance_ratio = n_accepted_mutants / n_generated_mutants

    mutation_number = [i for i in 0:n_accepted_mutants]
    stripe_indicator =  map(ft->ft[1] == 0,fitness_traj_tuple)
    H0 = minimum(findall(stripe_indicator))

    fitness_traj_add = map(ft->add_fitness(ft),fitness_traj_tuple)
    top_edits = compute_cumulative_edits(reduce(hcat,topologies))
    weight_edits = nothing 
    masked_hamming_distance_H0 = nothing
    masked_hamming_distance = nothing

    initial_fitness = fitness_traj_add[1]
    final_fitness = fitness_traj_add[end]

    ######## MST

    minimal_stripe_subgraphs = nothing 
    metagraph_vertices = nothing
    metagraph_parents = nothing

    ######## generate mutation data

    weight_id = map(mc-> [mut_choice_map[i] for i in mc],mut_choices)

    fitness_delta = get_fitness_delta(fitness_traj_add)

    n_weight_changes = map(m->length(m),weight_id)

    weight_id_label = map(v->join(map(x->weight_names[x],v),"|"),weight_id)

    start_fitness = fitness_traj_add[1:end-1]
    start_fitness_tuple = fitness_traj_tuple[1:end-1]
    mutant_fitness = fitness_traj_add[2:end]

    start_network = [collect(v) for v in eachcol(geno_traj_m[:,1:end-1])]

    is_new_int = [[sn[wi] == 0 for wi in mw_id] for (sn,mw_id) in zip(start_network,weight_id)]

    mut_types_combi = [[t for t in zip(new_int,mt)] for (new_int,mt) in zip(is_new_int,mut_types)]

    mutant_info = [(weight_id = weight_id[n],weight_id_label = weight_id_label[n],mut_type = mut_types_combi[n], new_interaction = is_new_int[n], mut_size = mut_sizes[n],start_fitness = start_fitness[n],start_fitness_tuple = start_fitness_tuple[n],mutant_fitness = mutant_fitness[n],fitness_delta = fitness_delta[n],start_network = start_network[n]) for n in 1:length(weight_id)]

    ####### predictions

    parent_inclusion_indicator = nothing

    tt_label_probabilities = nothing
    tt_label_predictions = nothing
    tt_label_entropies = nothing
    tt_prediction_error = nothing

    gt_label_probabilities = nothing
    gt_label_predictions = nothing
    gt_label_entropies = nothing
    gt_prediction_error = nothing

    mss_probabilities = nothing
    mss_predictions = nothing
    mss_entropies = nothing
    mss_prediction_error = nothing

    train_test_indicator = nothing
    train_test_indicator_mss = nothing
    embeddings = nothing

    epistasis = nothing

    tt_kl_div = nothing
    gt_kl_div = nothing
    mss_kl_div = nothing

    tt_shap = nothing
    gt_shap = nothing
    mss_shap = nothing

    other = nothing

    ####### instantiate 

    Trajectory(sim_id,geno_traj,topologies,n_accepted_mutants,acceptance_ratio,mutation_number,stripe_indicator,H0,wait_times_v,fitness_traj_add,fitness_traj_tuple,top_edits,weight_edits,masked_hamming_distance_H0,masked_hamming_distance,initial_fitness,final_fitness,mutant_info,metagraph_vertices,metagraph_parents,minimal_stripe_subgraphs,parent_inclusion_indicator,
                                                                                                            tt_label_probabilities,tt_label_predictions,tt_label_entropies,tt_prediction_error,gt_label_probabilities,gt_label_predictions,gt_label_entropies,gt_prediction_error,
                                                                                                            mss_probabilities,mss_predictions,mss_entropies,mss_prediction_error,train_test_indicator,train_test_indicator_mss,embeddings,epistasis,tt_kl_div,gt_kl_div,mss_kl_div,tt_shap,gt_shap,mss_shap,other)
end

function compute_cumulative_edits(gt)

    dham = pairwise(Hamming(),gt)

    total_ham = 0

    cumulative_ham = []

    push!(cumulative_ham,total_ham)

    for i in 1:size(gt,2)-1

        total_ham += dham[i,i+1]

        push!(cumulative_ham,total_ham)
    end

    return Int64.(cumulative_ham)

end

function masked_hamming_distance(topology,target_topology)

    new_topology = copy(topology)

    z0 = findall(x->x == 0,target_topology)

    new_topology[z0] .= 0.

    Distances.evaluate(Hamming(),new_topology,target_topology)

end

function assign_minimal_subgraphs!(tr::Trajectory,fs,ls)
    tr.minimal_stripe_subgraphs = fill([],length(tr.topologies))
    tr.minimal_stripe_subgraphs[tr.H0] = fs
    tr.minimal_stripe_subgraphs[end] = ls

    tr.masked_hamming_distance_H0 = [[masked_hamming_distance(top,Int.(sign.(fs))) for top in tr.topologies],[masked_hamming_distance(tr.topologies[1],top) for top in tr.topologies]]
    tr.masked_hamming_distance = [[masked_hamming_distance(top,Int.(sign.(ls))) for top in tr.topologies],[masked_hamming_distance(tr.topologies[1],top) for top in tr.topologies]]
end

function create_inclusion_metagraph(trajectories::Vector{Trajectory})

    min_stripe_top = unique(reduce(vcat,[tr.minimal_stripe_subgraphs[[tr.H0,end]] for tr in trajectories]))

    n_min_stripe_top = length(min_stripe_top)

    min_stripe_top_complexity = map(top->sum(abs.(top)),min_stripe_top)

    min_stripe_top_ordered = min_stripe_top[sortperm(min_stripe_top_complexity)]

    vertex_top_map = Dict(n=>top for (n,top) in enumerate(min_stripe_top_ordered));
    vertex_complexity_map = Dict(n=>sum(abs.(top)) for (n,top) in enumerate(min_stripe_top_ordered));
    top_vertex_map = Dict(top=>n for (n,top) in enumerate(min_stripe_top_ordered));

    inclusion_matrix = zeros(Int64,(n_min_stripe_top,n_min_stripe_top))

    for i in 1:n_min_stripe_top
        for j in 1:n_min_stripe_top
            if i != j
                inclusion_matrix[i,j] = test_inclusion(min_stripe_top_ordered[j],min_stripe_top_ordered[i])
            end
        end
    end

    inc_metagraph = SimpleDiGraph(inclusion_matrix)

    return inc_metagraph, vertex_top_map, top_vertex_map, vertex_complexity_map,inclusion_matrix
end


function assign_inc_vertex_ids!(tr::Trajectory,top_vertex_map)
    tr.inc_metagraph_vertices = [n ∈ [tr.H0,length(tr.topologies)] ? top_vertex_map[tr.minimal_stripe_subgraphs[n]] : -1 for n in 1:length(tr.topologies)]
end

function assign_inc_parents!(tr::Trajectory,inclusion_matrix,vertex_complexity_map,minimal_motif_id)

    tr.inc_metagraph_parents = []

    for n in 1:length(tr.topologies)

        if n ∈ [tr.H0,length(tr.topologies)]

            if tr.inc_metagraph_vertices[n] ∈ minimal_motif_id

                push!(tr.inc_metagraph_parents,tr.inc_metagraph_vertices[n])
            else
                options = minimal_motif_id[findall(inclusion_matrix[minimal_motif_id,tr.inc_metagraph_vertices[n]] .== 1)]
                choice = argmax([vertex_complexity_map[v] for v in options])

                push!(tr.inc_metagraph_parents,options[choice])
            end
        else
            push!(tr.inc_metagraph_parents,-1)
        end 
    end
end

function assign_weight_edits!(tr)

    we = vcat([0],map(mi->length(mi[:weight_id]),tr.mutant_info))

    c_we = [sum(we[1:i+1]) for i in 0:length(we)-1]

    tr.weight_edits = c_we
end

function select_minimal_topologies(list_mss)
    if length(list_mss) == 1
        return sign.(list_mss[1][3])
    else
        id = argmax(map(x->x[2],list_mss))
        return sign.(list_mss[id][3])
    end
end

const powerset_topologies = [[bit == '1' ? -1 : bit == '2' ? 0 : 1 for bit in string(i;base = 3,pad = 10)] for i in 0:3^10-1];


function v_restricted_cumulative_wait_time(tr::Trajectory,restriction,restriction_measure)
    id = maximum(findall(restriction_measure(tr,restriction)))

    return sum(tr.wait_times[1:id])
end

function cumulative_wait_time(tr::Trajectory)

    # return [sum(tr.wait_times[1:tr.H0-1])/(tr.H0-1);sum(tr.wait_times[tr.H0]);sum(tr.wait_times[tr.H0+1:end])/(length(tr.wait_times)-tr.H0-1)]

    a = tr.wait_times[1:tr.H0-1]
    b = tr.wait_times[tr.H0]
    c = tr.wait_times[tr.H0+1:end]

    if length(c) == 0 
        return [sum(a)/length(a);sum(b)/length(b);0]
    else
        return [sum(a)/length(a);sum(b)/length(b);sum(c)/length(c)]
    end
    
end

function average_wait_time(tr::Trajectory)

    # return [sum(tr.wait_times[1:tr.H0-1])/(tr.H0-1);sum(tr.wait_times[tr.H0]);sum(tr.wait_times[tr.H0+1:end])/(length(tr.wait_times)-tr.H0-1)]

    a = tr.wait_times[1:tr.H0-1]
    b = tr.wait_times[tr.H0]
    c = tr.wait_times[tr.H0+1:end]

    if length(c) == 0 
        return [sum(a);sum(b);0]
    else
        return [sum(a);sum(b);sum(c)]
    end
    
end

function average_wait_time_ext(tr::Trajectory)

    # return [sum(tr.wait_times[1:tr.H0-1])/(tr.H0-1);sum(tr.wait_times[tr.H0]);sum(tr.wait_times[tr.H0+1:end])/(length(tr.wait_times)-tr.H0-1)]

    a = tr.wait_times[2]
    b = tr.wait_times[2:tr.H0-1]
    c = tr.wait_times[tr.H0]
    d = tr.wait_times[tr.H0+1:end]

    all_v = []
    
    for v in [a,b,c,d]
        if length(v) == 0
            push!(all_v,0.)
        else
            push!(all_v,sum(v))
        end
    end

    return all_v
    
end

function fitness_restriction_measure(tr,restriction)
    tr.fitness_traj  .<= restriction
end

function calculate_probability(prob_time_slice,target_top,exact = true)

    p = 1.

    for (wn,w) in enumerate(target_top)
        if w == -1
            p *= prob_time_slice[1,wn]
        elseif w== 1
            p *= prob_time_slice[3,wn]
        else
            p *= prob_time_slice[2,wn]
        end
    end

    return p

end

function assign_predictions!(tr::Trajectory,model,prediction_type)

    if prediction_type == :tt
        tt = reduce(hcat,tr.topologies)
        tt_dtrain = xgboost.DMatrix(tt[1:10,:] |> transpose |> collect, feature_types = c_types, feature_names = weight_names)
        tr.tt_label_probabilities = model.predict(tt_dtrain)
        tr.tt_label_predictions = mapslices(p->argmax(p),tr.tt_label_probabilities,dims = 2)
        tr.tt_label_entropies = mapslices(p->entropy(p),tr.tt_label_probabilities,dims = 2);

    elseif prediction_type == :gt
        gt = reduce(hcat,tr.geno_traj)
        gt_dtrain = xgboost.DMatrix(gt[1:10,:] |> transpose |> collect, feature_names = weight_names)
        tr.gt_label_probabilities = model.predict(gt_dtrain)
        tr.gt_label_predictions = mapslices(p->argmax(p),tr.gt_label_probabilities,dims = 2)
        tr.gt_label_entropies = mapslices(p->entropy(p),tr.gt_label_probabilities,dims = 2);
    else
        gt = reduce(hcat,tr.geno_traj)
        gt_dtrain = xgboost.DMatrix(gt[1:10,:] |> transpose |> collect, feature_names = weight_names)

        prediction_prob = []
        prediction_labels = []

        for m in model
            edge_prob = m.predict(gt_dtrain)
            edge_label = mapslices(x->argmax(x)-2 ,edge_prob,dims = 2)

            push!(prediction_prob,edge_prob)
            push!(prediction_labels,edge_label)
        end

        tr.mss_probabilities = reduce((x,y) -> cat(x,y,dims = 3),[reshape(mss_p,(size(mss_p)...,1)) for mss_p in prediction_prob])
        tr.mss_predictions = [r |> collect for r in eachrow(reduce(hcat,prediction_labels))]
        
        mss_entropies = []

        for n in 1:length(tr.topologies)
            ps = @view tr.mss_probabilities[n,:,:]
            e = entropy([calculate_probability(ps,t) for t in powerset_topologies])
            push!(mss_entropies,e)
        end

        tr.mss_entropies = mss_entropies
    end

end

function assign_predictions!(tr::Trajectory,model,prediction_type,predict_label_to_vertex)

    if prediction_type == :tt
        tt = reduce(hcat,tr.topologies)
        tt_dtrain = xgboost.DMatrix(tt[1:10,:] |> transpose |> collect, feature_types = c_types, feature_names = weight_names)
        tr.tt_label_probabilities = model.predict(tt_dtrain)
        tr.tt_label_predictions = mapslices(p->predict_label_to_vertex[argmax(p)],tr.tt_label_probabilities,dims = 2)
        tr.tt_label_entropies = mapslices(p->entropy(p),tr.tt_label_probabilities,dims = 2);

    elseif prediction_type == :gt
        gt = reduce(hcat,tr.geno_traj)
        gt_dtrain = xgboost.DMatrix(gt[1:10,:] |> transpose |> collect, feature_names = weight_names)
        tr.gt_label_probabilities = model.predict(gt_dtrain)
        tr.gt_label_predictions = mapslices(p->predict_label_to_vertex[argmax(p)],tr.gt_label_probabilities,dims = 2)
        tr.gt_label_entropies = mapslices(p->entropy(p),tr.gt_label_probabilities,dims = 2);
    else
        gt = reduce(hcat,tr.geno_traj)
        gt_dtrain = xgboost.DMatrix(gt[1:10,:] |> transpose |> collect, feature_names = weight_names)

        prediction_prob = []
        prediction_labels = []

        for m in model
            edge_prob = m.predict(gt_dtrain)
            edge_label = mapslices(x->argmax(x)-2 ,edge_prob,dims = 2)

            push!(prediction_prob,edge_prob)
            push!(prediction_labels,edge_label)
        end

        tr.mss_probabilities = reduce((x,y) -> cat(x,y,dims = 3),[reshape(mss_p,(size(mss_p)...,1)) for mss_p in prediction_prob])
        tr.mss_predictions = [r |> collect for r in eachrow(reduce(hcat,prediction_labels))]
        
        mss_entropies = []

        for n in 1:length(tr.topologies)
            ps = @view tr.mss_probabilities[n,:,:]
            e = entropy([calculate_probability(ps,t) for t in powerset_topologies])
            push!(mss_entropies,e)
        end

        tr.mss_entropies = mss_entropies
    end

end

function assign_predictions_joint!(tr,model_label,all_models_mss,predict_label_to_vertex,top_vertex_map)

        gt = reduce(hcat,tr.geno_traj)
        gt_dtrain = xgboost.DMatrix(gt[1:10,:] |> transpose |> collect, feature_names = weight_names)
        tr.gt_label_probabilities = model_label.predict(gt_dtrain)
        tr.gt_label_predictions = mapslices(p->predict_label_to_vertex[argmax(p)],tr.gt_label_probabilities,dims = 2)
        tr.gt_label_entropies = mapslices(p->entropy(p),tr.gt_label_probabilities,dims = 2);

        gt_mss = vcat(gt[1:10,:],mapslices(p->argmax(p),tr.gt_label_probabilities,dims = 2) |> transpose |> collect)

        gt_mss_dtrain = xgboost.DMatrix(gt_mss[1:11,:] |> transpose |> collect, feature_names = weight_names_mss,feature_types=ft_mss, enable_categorical=true)

        prediction_prob = []
        prediction_labels = []

        for m in all_models_mss
            edge_prob = m.predict(gt_mss_dtrain)
            edge_label = mapslices(x->argmax(x)-2 ,edge_prob,dims = 2)

            push!(prediction_prob,edge_prob)
            push!(prediction_labels,edge_label)
        end

        tr.mss_probabilities = reduce((x,y) -> cat(x,y,dims = 3),[reshape(mss_p,(size(mss_p)...,1)) for mss_p in prediction_prob])
        tr.mss_predictions = [r |> collect for r in eachrow(reduce(hcat,prediction_labels))]
        
        mss_entropies = []

        for n in 1:length(tr.topologies)
            ps = @view tr.mss_probabilities[n,:,:]
            e = entropy([calculate_probability(ps,t) for t in powerset_topologies])
            push!(mss_entropies,e)
        end

        tr.mss_entropies = mss_entropies

end

function assign_predictions_sk!(tr::Trajectory,model,prediction_type,predict_label_to_vertex)

    if prediction_type == :tt
        tt = reduce(hcat,tr.topologies)
        tt_dtrain = tt[1:10,:] |> transpose |> collect
        tr.tt_label_probabilities = model.predict_proba(tt_dtrain)
        tr.tt_label_predictions = mapslices(p->predict_label_to_vertex[argmax(p)],tr.tt_label_probabilities,dims = 2)
        tr.tt_label_entropies = mapslices(p->entropy(p),tr.tt_label_probabilities,dims = 2);

    elseif prediction_type == :gt
        gt = reduce(hcat,tr.geno_traj)
        gt_dtrain = gt[1:10,:] |> transpose |> collect
        tr.gt_label_probabilities = model.predict_proba(gt_dtrain)
        tr.gt_label_predictions = mapslices(p->predict_label_to_vertex[argmax(p)],tr.gt_label_probabilities,dims = 2)
        tr.gt_label_entropies = mapslices(p->entropy(p),tr.gt_label_probabilities,dims = 2);
    else
        nothing
    end

end


function assign_predictions_sk!(tr::Trajectory,model,prediction_type,predict_label_to_vertex,weight_class_dict)

    if prediction_type == :tt
        tt = reduce(hcat,tr.topologies)
        tt_dtrain = tt[1:10,:] |> transpose |> collect
        tr.tt_label_probabilities = model.predict_proba(tt_dtrain)
        tr.tt_label_predictions = mapslices(p->predict_label_to_vertex[argmax(p)],tr.tt_label_probabilities,dims = 2)
        tr.tt_label_entropies = mapslices(p->entropy(p),tr.tt_label_probabilities,dims = 2);

    elseif prediction_type == :gt
        gt = reduce(hcat,tr.geno_traj)
        gt_dtrain = gt[1:10,:] |> transpose |> collect
        tr.gt_label_probabilities = model.predict_proba(gt_dtrain)
        tr.gt_label_predictions = mapslices(p->predict_label_to_vertex[argmax(p)],tr.gt_label_probabilities,dims = 2)
        tr.gt_label_entropies = mapslices(p->entropy(p),tr.gt_label_probabilities,dims = 2);
    else
        gt = reduce(hcat,tr.geno_traj)
        gt_dtrain = gt[1:10,:] |> transpose |> collect

        prediction_prob = []
        prediction_labels = []

        for (n,m) in enumerate(model)
            edge_prob_v = m.predict_proba(gt_dtrain)

            edge_prob = zeros(size(edge_prob_v,1),3)

            for (nc,col) in enumerate(eachcol(edge_prob_v))
                edge_prob[:,weight_class_dict[n][2][nc-1]+2] = col
            end

            edge_label = mapslices(x->argmax(x)-2 ,edge_prob,dims = 2)

            push!(prediction_prob,edge_prob)
            push!(prediction_labels,edge_label)
        end

        tr.mss_probabilities = reduce((x,y) -> cat(x,y,dims = 3),[reshape(mss_p,(size(mss_p)...,1)) for mss_p in prediction_prob])
        tr.mss_predictions = [r |> collect for r in eachrow(reduce(hcat,prediction_labels))]
        
        mss_entropies = []

        for n in 1:length(tr.topologies)
            ps = @view tr.mss_probabilities[n,:,:]
            e = entropy([calculate_probability(ps,t) for t in powerset_topologies])
            push!(mss_entropies,e)
        end

        tr.mss_entropies = mss_entropies
    end

end

function assign_mss_prediction_errors!(tr::Trajectory,metric)

    tr.mss_prediction_error = [Distances.evaluate(metric,vcat(pred,[0.,0.]),tr.minimal_stripe_subgraphs[tr.H0]) for pred in tr.mss_predictions]

end

function assign_mss_prediction_errors!(tr::Trajectory)

    tr.mss_prediction_error = [vcat(pred,[0.,0.]) == tr.minimal_stripe_subgraphs[tr.H0] for pred in tr.mss_predictions]

end

function assign_tt_prediction_errors!(tr::Trajectory,label)

    tr.tt_prediction_error = [pred == label for pred in tr.tt_label_predictions]

end

function assign_gt_prediction_errors!(tr::Trajectory,label)

    tr.gt_prediction_error = [pred == label for pred in tr.gt_label_predictions]

end

function assign_tt_other_prediction_errors!(tr::Trajectory,label,predict_id)

    tr.tt_prediction_error = [label == -1 ? !(pred ∈ predict_id) :  pred == label for pred in tr.tt_label_predictions]

end

function assign_gt_other_prediction_errors!(tr::Trajectory,label,predict_id)

    tr.gt_prediction_error = [label == -1 ? !(pred ∈ predict_id) :  pred == label for pred in tr.gt_label_predictions]

end

function v_restricted_accuracy(tr::Trajectory,restriction_measure,prediction_type)
    prediction_id = maximum(findall(restriction_measure(tr)))
    if prediction_type == :tt
        return tr.tt_prediction_error[prediction_id]
    elseif prediction_type == :gt
        return tr.gt_prediction_error[prediction_id]
    else
        return tr.mss_prediction_error[prediction_id]
    end
end

function v_restricted_entropy(tr::Trajectory,restriction_measure,prediction_type)
    prediction_id = maximum(findall(restriction_measure(tr)))
    if prediction_type == :tt
        return tr.tt_label_entropies[prediction_id]
    elseif prediction_type == :gt
        return tr.gt_label_entropies[prediction_id]
    else
        return tr.mss_entropies[prediction_id]
    end
end

function v_restricted_probabilities(tr::Trajectory,restriction_measure,prediction_type)
    prediction_id = maximum(findall(restriction_measure(tr)))
    if prediction_type == :tt
        return tr.tt_label_probabilities[prediction_id,:]
    elseif prediction_type == :gt
        return tr.gt_label_probabilities[prediction_id,:]
    else
        return tr.mss_probabilities[prediction_id,:]
    end
end

function v_restricted_label_inclusion(tr::Trajectory,restriction_measure,label_type)
    prediction_id = maximum(findall(restriction_measure(tr)))
    if label_type == :H0
        return tr.H0 <= prediction_id
    else
        return length(tr.topologies) <= prediction_id
    end
end

function top_edit_restriction_measure(tr,restriction)
    tr.top_edits .<= restriction
end

function weight_edit_restriction_measure(tr,restriction)
    tr.weight_edits.<= restriction
end

function relative_fitness_restriction_measure(tr,restriction,label_type)
    if label_type == :H0
        ((tr.fitness_traj .-  tr.fitness_traj[1]) ./ tr.fitness_traj[tr.H0-1]) .<= restriction
    else
        ((tr.fitness_traj .-  tr.fitness_traj[1]) ./ tr.fitness_traj[end-1]) .<= restriction
    end
end

function relative_top_edit_restriction_measure(tr,restriction,label_type)
    if label_type == :H0
        (tr.top_edits ./ tr.top_edits[tr.H0]) .<= restriction
    else
        (tr.top_edits./ tr.top_edits[end]) .<= restriction
    end
end

function relative_top_edit_restriction_measure(tr,restriction,label_type)
    if label_type == :H0
        (tr.top_edits ./ tr.top_edits[tr.H0]) .<= restriction
    else
        (tr.top_edits./ tr.top_edits[end]) .<= restriction
    end
end

function masked_hamming_restriction_measure(tr,restriction,label_type)
    if label_type == :H0
        tr.masked_hamming_distance_H0[1] .>= restriction
    else
        tr.masked_hamming_distance[1] .>= restriction
    end
end

function masked_hamming_restriction_measure_begin(tr,restriction)

    tr.masked_hamming_distance[2] .<= restriction

end

function filter_by_value_and_train(tr,value,train_type)

    valid = findall(tr.top_edits .== value) ∩ findall(tr.train_test_indicator .== train_type)

    if length(valid) > 0

        return true
    else
        return false
    end

end

function filter_by_value(tr,value)

    valid = findall(tr.top_edits .== value)

    if length(valid) > 0

        return true
    else
        return false
    end

end