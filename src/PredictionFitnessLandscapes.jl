function set_weight(entry::Tuple{Int,Int},step::Float64,w::Matrix{Float64})
    new_w = copy(w)
    new_w[entry...] = step
    return new_w
end

function set_weight(w,entry_1::Tuple{Int,Int},entry_2::Tuple{Int,Int},weight_1::Float64,weight_2::Float64)
    new_w = copy(w)

    new_w[entry_1...] = weight_1
    new_w[entry_2...] = weight_2

    return vec(new_w)
end

function create_mutant_get_fitness(founder::Individual,development::DESystemSolver,entry_1::Tuple{Int,Int},entry_2::Tuple{Int,Int},weight_1::Float64,weight_2::Float64,fitness_function)

    new_w = copy(founder.genotype.p[1])

    new_w[entry_1...] = weight_1
    new_w[entry_2...] = weight_2

    mutant = Individual(remake(founder.genotype, p = (new_w,founder.genotype.p[2:end]...)),development)

    mutant_fitness = fitness_function(mutant.phenotype)

    return mutant_fitness, vec(new_w)

end

function assign_mst_label(v,top_vertex_map,predict_id,vertex_to_predict_label)
    if isnan(v[1])
        return 0
    else
        if haskey(top_vertex_map,v)
            if top_vertex_map[v] ∈ predict_id
                return vertex_to_predict_label[top_vertex_map[v]]
            else
                return 5
            end
        else
            return 5
        end
    end
end

function get_max_prob_point(n1,n2,pred_grid,prob_grid,fitness_grid,pred_choice,fitness_tol)

    fitness_valid = findall(x->x>fitness_tol,fitness_grid[n1,n2,:,:])
    pred_valid = findall(x->x==pred_choice,pred_grid[n1,n2,:,:])

    valid = fitness_valid ∩ pred_valid

    if length(valid) > 0
        prob_id = argmax(prob_grid[n1,n2,:,:][valid])

        return valid[prob_id]
    else
        return NaN
    end

end

function get_max_prob_point(n1,n2,pred_grid,prob_grid,fitness_grid,pred_choice,fitness_tol,top_choice,sample_grid_v)

    fitness_valid = first.(Tuple.(findall(x->x>fitness_tol,fitness_grid[n1,n2,:,:])))
    pred_valid = first.(Tuple.(findall(x->x==pred_choice,pred_grid[n1,n2,:,:])))
    top_valid_1 = findall(x->sign(x) == sign(top_choice[1]),sample_grid_v[1,:])
    top_valid_2 = findall(x->sign(x) == sign(top_choice[2]),sample_grid_v[2,:])

    valid = fitness_valid ∩ pred_valid ∩ top_valid_1 ∩ top_valid_2

    if length(valid) > 0
        prob_id = argmax(prob_grid[n1,n2,:,:][valid])

        return valid[prob_id]
    else
        print(n1)
        print(n2)
        return NaN
    end

end
