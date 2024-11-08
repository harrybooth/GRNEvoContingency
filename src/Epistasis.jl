function characterise_epistasis_deltas(combi_result,fitness_eps) # main

    if length(combi_result) > 2
        ratio_new_mutant = combi_result ./ (combi_result[end] - fitness_eps)

        accept_new_mutant = ratio_new_mutant[2:end-1] .>= 1

        if !any(accept_new_mutant)
            rtype = :rse
        elseif all(accept_new_mutant)
            rtype = :ne
        else
            rtype = :se
        end

        return rtype,combi_result
    else
        return :sm, combi_result
    end
end

function evaluate_epistasis_class_full_deltas(mut_tuple,grn_parameters,development,fitness_function,mut_op::MutationOperatorDual) # main

    n_mut = length(mut_tuple[:weight_id])

    if n_mut > 1

        mut_combi = [[bit == '1' ? true : false for bit in string(i;base = 2,pad = n_mut)] for i in 1:2^n_mut-1]

        prob_mutant = []

        for mut_id in mut_combi[1:end]

            new_network = noise_specified(mut_tuple[:start_network],mut_tuple[:weight_id][mut_id],mut_tuple[:mut_size][mut_id],mut_tuple[:mut_type][mut_id],mut_op)

            mutant = Individual(reshape(new_network,(3,4)),grn_parameters,development)

            mutant_fitness = fitness_function(mutant.phenotype)

            Δf1 = mutant_fitness[1] - mut_tuple[:start_fitness_tuple][1]
            Δf2 = mutant_fitness[2] - mut_tuple[:start_fitness_tuple][2]

            Δf = Δf1 != 0. ? Δf1 : Δf2

            push!(prob_mutant,Δf)

        end

        return vcat([0.],prob_mutant)

    else
        return [0.,1.]
    end
end

function calculate_epi_class_proportion(epi_class)
    
    total = length(epi_class)

    return [count(x->x==ec,epi_class)/total for ec in [:rse,:se,:ne,:sm]]
end

function calculate_epi_class_proportion_fw(epi_class,fitness_deltas)
    
    total_fitness_change = sum(fitness_deltas)

    return [sum(fitness_deltas[findall(x->x==ec,epi_class)])/total_fitness_change for ec in [:rse,:se,:ne,:sm]]
end

function calculate_epi_class_proportion(epi_class,epi_class_totals)
    
    return [count(x->x==ec,epi_class) for ec in [:rse,:se,:ne,:sm]] ./ epi_class_totals
end

function evaluate_epistasis_types!(tr,grn_parameters,development,fitness_function)
    all_class_epi = map(mi->evaluate_epistasis_class(mi,grn_parameters,development,fitness_function),tr.mutant_info);
    tr.epistasis = all_class_epi
end

function evaluate_epistasis_types!(tr,grn_parameters,development,fitness_function,mut_op::MutationOperatorDual)
    all_class_epi = map(mi->evaluate_epistasis_class(mi,grn_parameters,development,fitness_function,mut_op),tr.mutant_info);
    tr.epistasis = all_class_epi
end

function evaluate_epistasis_types_shap!(tr,grn_parameters,development,fitness_function,mut_op::MutationOperatorDual,fitness_eps)
    all_class_epi = map(mi->evaluate_epistasis_class_shapley(mi,grn_parameters,development,fitness_function,mut_op,fitness_eps),tr.mutant_info);
    tr.epistasis = all_class_epi
end

function evaluate_epistasis_types_full!(tr,grn_parameters,development,fitness_function,mut_op::MutationOperatorDual)
    all_class_epi = map(mi->evaluate_epistasis_class_full(mi,grn_parameters,development,fitness_function,mut_op),tr.mutant_info);
    tr.epistasis = all_class_epi
end

function evaluate_epistasis_types_full_deltas!(tr,grn_parameters,development,fitness_function,mut_op::MutationOperatorDual)
    all_class_epi = map(mi->evaluate_epistasis_class_full_deltas(mi,grn_parameters,development,fitness_function,mut_op),tr.mutant_info);
    tr.epistasis = all_class_epi
end

function get_mut_size_by_type(tr,type,range_l,range_u)

    all_type = []

    for mi in tr.mutant_info[range_l:range_u]

        mi_id = findall(x->x==type,mi.mut_type)

        if length(mi_id) >= 1
            push!(all_type,mi[:mut_size][mi_id])
        end

    end

    if length(all_type) > 0

        return reduce(vcat,all_type)
    else
        return []
    end

end

function get_mut_size_by_type_and_weight(tr,type,weight_id,range_l,range_u)

    all_type = []

    for mi in tr.mutant_info[range_l:range_u]

        mi_id = findall(x->x==type,mi.mut_type)
        mi_weight_id = findall(x->x==weight_id,mi.weight_id)

        look_id = mi_id ∩ mi_weight_id

        if length(look_id) >= 1
            push!(all_type,mi[:mut_size][look_id])
        end
    end

    if length(all_type) > 0

        return reduce(vcat,all_type)
    else
        return []
    end

end

function get_mut_type(tr,range_l,range_u)

    d = [mi[:mut_type] for mi in tr.mutant_info[range_l:range_u]]

    if length(d) >= 1
        return reduce(vcat,d)
    else
        return []
    end

end

function get_mut_type_full(tr,range_l,range_u)

    d = [mi[:mut_type] for mi in tr.mutant_info[range_l:range_u]]

    if length(d) >= 1
        return d
    else
        return []
    end

end

function get_mut_n(tr,range_l,range_u)

    d = [length(mi[:weight_id]) for mi in tr.mutant_info[range_l:range_u]]

    if length(d) >= 1
        return reduce(vcat,d)
    else
        return []
    end

end

function calculate_mut_type_proportion(mut_type,types)
    
    total = length(mut_type)

    return [count(x->x==ec,mut_type)/total for ec in types]
end

function calculate_mut_type_proportion_list(mut_type,types)
    
    total = length(mut_type)

    return [count(x->x ∈ ec,mut_type)/total for ec in types]
end

function calculate_mut_type_count(mut_type,types)

    return [count(x->x==ec,mut_type) for ec in types]
end

function nothing_new_in_common(target_top,new_network,initial_network)

    result = true

    for i in 1:10
        if (initial_network[i] == 0)
            if (target_top[i] != 0) && (sign(target_top[i]) == sign(new_network[i])) 
                result = false
            end
        end
    end

    return result

end 

function get_mut_weight_full(tr,range_l,range_u,weight_loop_dict)

    d = [[weight_loop_dict[x] for x in mi[:weight_id]] for mi in tr.mutant_info[range_l:range_u]]

    if length(d) >= 1
        return d
    else
        return []
    end

end

function get_mut_weight_full_new(tr,range_l,range_u,weight_loop_dict)

    d = [[(weight_loop_dict[x],t[1]) for (x,t) in zip(mi[:weight_id],mi[:mut_type])] for mi in tr.mutant_info[range_l:range_u]]

    if length(d) >= 1
        return d
    else
        return []
    end

end
