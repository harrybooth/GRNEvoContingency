
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

function get_rel_dyn_vector(n1,t1,n_steps,save_id)

    w_ind = Individual(n1,t1,grn_parameters,TimeStampedGRNSolver(LinRange(0,t1,n_steps),save_id))

    return reduce(vcat,w_ind.phenotype.u[2:end-1])

end

function get_av_dyn_vector(n1,t1,n_steps,n_segments)

    w_ind = Individual(n1,t1,grn_parameters,TimeStampedGRNSolver(LinRange(0,t1,n_steps)))

    return reduce(vcat,map(x->vec(reduce(hcat,[mean(x[:,t:t+Int(Nc/n_segments)-1],dims = 2) for t in 1:Int(Nc/n_segments):Nc])),w_ind.phenotype.u[2:end-1]))

end

function get_dense_dyn_sol(n1,t1,save_id)

    w_ind = Individual(n1,t1,grn_parameters,DenseGRNSolver(save_id));

    ind1_f = t->w_ind.phenotype(t*t1)

    return ind1_f

end

##### Assignment tools

function test_inclusion(net_v,top_v)

    n1 = sign.(net_v)
    incl = 1

    for (n,w) in enumerate(top_v)

        if w != 0

            if n1[n] != w
                incl = 0
            end
        end
    end

    return incl
end

function assign_class(row)
    if all(row.==0)
        return "no assignment"
    else
        return fundamental_topologies[findall(x->x==1,row)]
    end
end

function determine_class(en_top,dyn_top)

    r = zeros(Int,size(en_top,1))

    for i in 1:size(en_top,1)
        id =  findall(x->x==1,en_top[i,:])
        if length(id) == 0
            r[i] = 0
        elseif length(id) == 1
            r[i] = id[1]
        else
            comb = en_top[i,:] .& dyn_top[i,:]

            if all(comb.==0)
                r[i] = 0
            else
                r[i] = findall(x->x==1,comb)[1]
            end
        end
    end

    return r
end

function entropy(v)
    length(v) != 1 ? -sum(v[v .!= 0] .* log.(v[v .!= 0])) / log(length(v)) : -sum(v[v .!= 0] .* log.(v[v .!= 0]))
end

function cross_entropy(p,q)

    total_entropy = 0

    for (pi,qi) in zip(p,q)
        if (pi == 0) || (qi == 0)
            nothing
        else
            total_entropy += pi*log(qi)
        end
    end
    
    return -total_entropy
end
