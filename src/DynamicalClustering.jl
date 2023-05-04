using DynamicAxisWarping

function calculate_dyn_distance(n1,t1,n2,t2,dist,n_steps,save_id,relative)

    if relative
        w_ind_1 = Individual(n1,t1,grn_parameters,TimeStampedGRNSolver(LinRange(0,t1,n_steps),save_id));
        w_ind_2 = Individual(n2,t2,grn_parameters,TimeStampedGRNSolver(LinRange(0,t2,n_steps),save_id));
    else
        tf = min(t1,t2)
        w_ind_1 = Individual(n1,tf,grn_parameters,TimeStampedGRNSolver(LinRange(0,tf,n_steps),save_id));
        w_ind_2 = Individual(n2,tf,grn_parameters,TimeStampedGRNSolver(LinRange(0,tf,n_steps),save_id));
    end

    v1 = reduce(vcat,w_ind_1.phenotype.u[2:end-1])
    v2 = reduce(vcat,w_ind_2.phenotype.u[2:end-1])

    return evaluate(dist, v1, v2)
end

function calculate_dyn_distance(n1,t1,w_ind_2,dist,n_steps,save_id)

    w_ind_1 = Individual(n1,t1,grn_parameters,TimeStampedGRNSolver(LinRange(0,t1,n_steps),save_id));

    v1 = reduce(vcat,w_ind_1.phenotype.u[2:end-1])
    v2 = reduce(vcat,w_ind_2.phenotype.u[2:end-1])

    return evaluate(dist, v1, v2)
end

function calculate_dyn_dtw(n1,n2,t1,t2,save_id)

    w_ind_1 = Individual(n1,t1,grn_parameters,DenseGRNSolver(save_id));
    w_ind_2 = Individual(n2,t2,grn_parameters,DenseGRNSolver(save_id));

    ind1_f = t->w_ind_1.phenotype(t*t1)
    ind2_f = t->w_ind_2.phenotype(t*t2)

    cost, ϕ, ψ = gdtw(ind1_f,ind2_f)

    return cost
end

function get_rel_dyn_vector(n1,t1,n_steps,save_id)

    w_ind = Individual(n1,t1,grn_parameters,TimeStampedGRNSolver(LinRange(0,t1,n_steps),save_id))

    return reduce(vcat,w_ind.phenotype.u[2:end-1])

end

function get_dense_dyn_sol(n1,t1,save_id)

    w_ind = Individual(n1,t1,grn_parameters,DenseGRNSolver(save_id));

    ind1_f = t->w_ind.phenotype(t*t1)

    return ind1_f

end

function calculate_dyn_dtw(ind1_f,ind2_f)

    cost, ϕ, ψ = gdtw(ind1_f,ind2_f)

    return cost
end