function create_label_H0_parent(tr,top_vertex_map,predict_id,predict_mst,predict_mst_complexity)
    if tr.minimal_stripe_subgraphs[tr.H0] ∈ predict_mst
        # vertex_to_predict_label[tr.inc_metagraph_vertices[tr.H0]]
        top_vertex_map[tr.minimal_stripe_subgraphs[tr.H0]]
    else
        incl = [test_inclusion(tr.minimal_stripe_subgraphs[tr.H0],pmst) for pmst in predict_mst]

        if sum(incl) != 0
            inc_id = findall(x->x==1, incl)
            choice = inc_id[argmax(predict_mst_complexity[inc_id])]
            predict_id[choice]
        else
            -1
        end

    end
end

function create_label_H0(tr,top_vertex_map,predict_id,predict_mst)
    if tr.minimal_stripe_subgraphs[tr.H0] ∈ predict_mst
        top_vertex_map[tr.minimal_stripe_subgraphs[tr.H0]]
    else
        -1
    end
end