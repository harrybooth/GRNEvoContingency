function make_square(m::Vector)
    vcat(reshape(m,(3,4)),zeros(typeof(m[1]),(1,4)))
end

function add_fitness(tuple_f)
    return tuple_f[1] + ((tuple_f[2]+1)/2) + 1
end

function create_full_fitness_traj(fitness_traj,wait_times)

    all_ff = []
    
    for (fitness,wt) in zip(fitness_traj[1:end-1],wait_times)

        full_fitness  = fill(fitness,wt)

        push!(all_ff,full_fitness)

    end

    all_ff_v = reduce(vcat,all_ff)

    all_ff_vf = vcat(all_ff_v,[fitness_traj[end]])

    return all_ff_vf
end


function group_values_into_bins(values,bins)
    bin_counts = zeros(Int,length(bins)-1)

    new_values = copy(values)

    for i in 1:length(bins)-1
        total = 0.
        for (n,v) in enumerate(new_values)
            if (v < bins[i+1]) && (v >= bins[i])
                total+=1 
                deleteat!(new_values,n)
            end
        end
        bin_counts[i] = total
    end

    return bin_counts
end

uniqueidx(x) = unique(i -> x[:,i], 1:size(x,2))

uniqueid(x) = unique(i -> x[i], 1:length(x))
    

function cond_save(dir,fig,cond)
    if cond
        CairoMakie.save(dir,fig,pt_per_unit = 1)
    end
end

function return_order_by_count(v)

    v_un = unique(v)
    counts_v = [count(x->x==value,v) for value in v_un]

    order_v = sortperm(counts_v,rev = true)

    return v_un[order_v],counts_v[order_v]
end