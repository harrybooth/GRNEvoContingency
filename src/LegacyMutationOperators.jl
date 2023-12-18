struct MutationOperatorSeq 
    min_affinity :: Any
    max_affinity :: Any
    ϵ :: Any
    Ls :: Any
    μ :: Any
    flip_prob :: Any
    mutation_weights :: Vector{CartesianIndex{2}}
end

struct MutationOperatorSeqApprox
    noise_distribution :: Distribution
    n_sample_func :: Any
    start_affinity :: Any
    max_affinity :: Any
    flip_prob :: Any
    mutation_weights :: Vector{CartesianIndex{2}}
end

mutable struct TF
    name
    associated_w
    mutate_probability
    pm_probability
    pm_noise_distribution
    duplication_noise_distribution
    start_affinity 
    reg_flip_probability
end

mutable struct TFBS
    name
    associated_w
    mutate_probability
    pm_probability
    pm_noise_distribution
    duplication_noise_distribution
    start_affinity
    reg_flip_probability
end

function noise(w::Matrix{Float64},mut_op::MutationOperator)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    choices = sample(mut_op.mutation_weights,n_mut,replace = false)

    for index in choices
        if new_w[index] == 0
            # proposal = new_w[index] + rand(mut_op.noise_distribution)
            new_w[index] = rand(Uniform(-mut_op.max_w,mut_op.max_w))
        else
            if rand() < mut_op.deletion_p
                del_index = sample(mut_op.mutation_weights,1)[1]
                new_w[del_index] = 0.

                proposal = new_w[index] + rand(mut_op.noise_distribution)*new_w[index]
                new_w[index] = abs(proposal) > mut_op.max_w ? mut_op.max_w*sign(proposal) : proposal
            else
                proposal = new_w[index] + rand(mut_op.noise_distribution)*new_w[index]
                new_w[index] = abs(proposal) > mut_op.max_w ? mut_op.max_w*sign(proposal) : proposal
            end
        end
    end

    valid = maximum(abs.(new_w)) <= mut_op.max_w

    return new_w,nothing,nothing,nothing, valid
end

function noise_mtype_seq_approx(w::Matrix{Float64},mut_op::MutationOperatorSeqApprox)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    choices = sort(sample(mut_op.mutation_weights,n_mut,replace = false))
    mtype = []
    sizes = []

    for index in choices
        if new_w[index] == 0
            n = exp(-1*rand(mut_op.noise_distribution))
            push!(mtype,:new)

            if rand() < 0.5
                new_w[index] = -1*mut_op.start_affinity*n
                push!(sizes,-n)
            else
                new_w[index] = mut_op.start_affinity*n
                push!(sizes,n)
            end
        else
            n = exp(-1*rand(mut_op.noise_distribution))

            new_w[index] = new_w[index]*n

            push!(mtype,:existing)

            if rand() < mut_op.flip_prob
                new_w[index] = -1*new_w[index]
                push!(sizes,-n)
            else
                push!(sizes,n)
            end
        end
    end

    valid = maximum(abs.(new_w)) <= mut_op.max_affinity

    return new_w, choices, mtype, sizes, valid
end

function noise_mtype_sequence(w::Matrix{Float64},mut_op::MutationOperatorSeq)

    new_w = copy(w)

    mtype = []
    sizes = []

    for index in mut_op.mutation_weights

        current_weight = w[index]

        d0 = (log(mut_op.max_affinity)-log(current_w))/mut_op.ϵ

        Z = Skellam((mut_op.Ls-d0)*mut_op.μ, d0*mut_op.μ/(mut_op.k-1))

        if rand(Z) != 0

            # m_compound = ((mut_op.Ls-d0)*mut_op.μ - d0*mut_op.μ/(mut_op.k-1))
            # v_compound = ((mut_op.Ls-d0)*μ + d0*mut_op.μ/(mut_op.k-1))

            Z_approx = truncated(Normal(mean(Z),var(Z)); lower = -d0, upper = mut_op.Ls - d0)

            n =  exp.(-ϵ .* rand(Z_approx))

            if rand() < mut_op.flip_prob
                new_w[index] = -current_weight * n
                push!(sizes,-n)
            else
                new_w[index] = current_weight * n
                push!(sizes,n)
            end

            if (current_weight <= mut_op.min_affinity) & (new_w[index] > mut_op.min_affinity)
                push!(mtype,:new)
            elseif (current_weight > mut_op.min_affinity) & (new_w[index] <= mut_op.min_affinity)
                push!(mtype,:del)
            elseif (current_weight > mut_op.min_affinity) & (new_w[index] > mut_op.min_affinity)
                push!(mtype,:existing)
            else
                push!(mtype,:other)
            end
        end
    end

    valid = maximum(abs.(new_w)) <= mut_op.max_affinity

    return new_w, choices, mtype, sizes, valid
end

function noise_mtype_dup(w::Matrix{Float64},all_sites,n_sample_func)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = n_sample_func()
    end

    choices = sample(all_sites,n_mut,replace = false)
    mtype = []
    all_sizes = []

    for site in choices

        sizes = []

        if rand() < site.pm_probability

            if typeof(site) == TF
                push!(mtype,:TF_pm)
                for index in site.associated_w
                    n = exp(-1*rand(site.pm_noise_distribution))
    
                    if new_w[index] == 0
            
                        if rand() < 0.5
                            new_w[index] = -1*site.start_affinity*n
                            push!(sizes,-n)
                        else
                            new_w[index] = site.start_affinity*n
                            push!(sizes,n)
                        end
                    else
    
                        if rand() < site.reg_flip_probability
                            new_w[index] = -1*n*new_w[index]
                            push!(sizes,-n)
                        else
                            new_w[index] = n*new_w[index]
                            push!(sizes,n)
                        end
                    end
                end
            else
                push!(mtype,:TFBS_pm)
                n = exp(-1*rand(site.pm_noise_distribution))
    
                if new_w[site.associated_w] == 0
        
                    if rand() < 0.5
                        new_w[site.associated_w] = -1*site.start_affinity*n
                        push!(sizes,-n)
                    else
                        new_w[site.associated_w] = site.start_affinity*n
                        push!(sizes,n)
                    end
                else

                    if rand() < site.reg_flip_probability
                        new_w[site.associated_w] = -1*n*new_w[site.associated_w]
                        push!(sizes,-n)
                    else
                        new_w[site.associated_w] = n*new_w[site.associated_w]
                        push!(sizes,n)
                    end
                end
            end

        else
            if typeof(site) == TF
                push!(mtype,:TF_dup)
                if site.name == "A"
                    target_gene = sample([2,3],1)[1]
                elseif site.name == "B"
                    target_gene = sample([1,3],1)[1]
                else
                    target_gene = sample([1,2],1)[1]
                end

                for (ng,index) in enumerate(site.associated_w)
                    n = new_w[index]*exp(-1*rand(site.pm_noise_distribution))
                    new_w[ng,target_gene] = new_w[ng,target_gene] + n
                    push!(sizes,n)
                end

            else
                push!(mtype,:TFBS_dup)
                target_site = sample(all_sites[4:end],1,replace = false)[1]
                n = new_w[site.associated_w]*exp(-1*rand(site.pm_noise_distribution))
                new_w[target_site.associated_w] = new_w[target_site.associated_w] + n
                push!(sizes,n)
            end
        end

        push!(all_sizes,sizes)
    end

    valid = maximum(abs.(new_w)) <= 10.

    return new_w, map(x->x.name,choices), mtype, all_sizes, valid
end

function noise_mtype_dup_1(w::Matrix{Float64},all_sites,n_sample_func)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = n_sample_func()
    end

    choices = sample(all_sites,n_mut,replace = false)
    mtype = []
    all_sizes = []

    for site in choices

        sizes = []

        if rand() < site.pm_probability

            if typeof(site) == TF
                push!(mtype,:TF_pm)
                for index in site.associated_w
                    n = exp(-1*rand(site.pm_noise_distribution))
    
                    if new_w[index] == 0
                        nothing
                    else
                        if rand() < site.reg_flip_probability
                            new_w[index] = -1*n*new_w[index]
                            push!(sizes,-n)
                        else
                            new_w[index] = n*new_w[index]
                            push!(sizes,n)
                        end
                    end
                end
            else
                push!(mtype,:TFBS_pm)
                n = exp(-1*rand(site.pm_noise_distribution))
    
                if new_w[site.associated_w] == 0
                    nothing
                else
                    if rand() < site.reg_flip_probability
                        new_w[site.associated_w] = -1*n*new_w[site.associated_w]
                        push!(sizes,-n)
                    else
                        new_w[site.associated_w] = n*new_w[site.associated_w]
                        push!(sizes,n)
                    end
                end
            end

        else
            if typeof(site) == TF

                push!(mtype,:TF_dup)

                target_gene = sample([1,2,3],1)[1]
  
                for (ng,index) in enumerate(site.associated_w)

                    if rand() < site.reg_flip_probability
                        n = -1*new_w[index]*exp(-1*rand(site.pm_noise_distribution))
                    else
                        n = -new_w[index]*exp(-1*rand(site.pm_noise_distribution))
                    end

                    new_w[ng,target_gene] = new_w[ng,target_gene] + n
                    push!(sizes,n)
                end

            else
                if site.name == "M=>A"
                    nothing
                else
                    push!(mtype,:TFBS_dup)
                    
                    target_gene = sample([1,2,3],1)[1]

                    if rand() < site.reg_flip_probability
                        n = -1*new_w[site.associated_w]*exp(-1*rand(site.pm_noise_distribution))
                    else
                        n = new_w[site.associated_w]*exp(-1*rand(site.pm_noise_distribution))
                    end

                    new_w[target_gene,site.associated_w[2]]  = new_w[target_gene,site.associated_w[2]] + n
                    push!(sizes,n)
                end
            end
        end

        push!(all_sizes,sizes)
    end

    valid = maximum(abs.(new_w)) <= 10.

    return new_w, map(x->x.associated_w,choices), mtype, all_sizes, valid
end

function noise_mtype_fm(w::Matrix{Float64},mut_op::MutationOperator,index)

    new_w = copy(w)

    choices = [CartesianIndex(index)]
    mtype = []
    sizes = []

    for index in choices
        if new_w[index] == 0
            n = rand(Uniform(-mut_op.max_w,mut_op.max_w))
            new_w[index] = n
            push!(mtype,:new)
            push!(sizes,n)
        else
            if rand() < mut_op.deletion_p
                push!(sizes,- new_w[index])
                new_w[index] = 0.
                push!(mtype,:del)
            else
                n = rand(mut_op.noise_distribution)
                new_w[index] = new_w[index] + n*new_w[index]
                push!(mtype,:existing)
                push!(sizes,n)
            end
        end
    end

    valid = maximum(abs.(new_w)) <= mut_op.max_w

    return new_w, choices, mtype, sizes, valid
end

function noise_runiform(w::Matrix{Float64},mut_op::MutationOperator,p_uniform)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    choices = sample(mut_op.mutation_weights,n_mut,replace = false)

    for index in choices
        if rand() < p_uniform
            new_w[index] = rand(Uniform(-mut_op.max_w,mut_op.max_w))
        else
            if new_w[index] == 0
                proposal = new_w[index] + rand(mut_op.noise_distribution)
                new_w[index] = abs(proposal) > mut_op.max_w ? mut_op.max_w*sign(proposal) : proposal
            else
                if rand() < mut_op.deletion_p
                    new_w[index] = 0.
                else
                    proposal = new_w[index] + rand(mut_op.noise_distribution)*new_w[index]
                    new_w[index] = abs(proposal) > mut_op.max_w ? mut_op.max_w*sign(proposal) : proposal
                end
            end
        end
    end

    valid = maximum(abs.(new_w)) <= mut_op.max_w

    return new_w,nothing,nothing,nothing, valid
end

function noise_no_additions(w::Matrix{Float64},mut_op::MutationOperator)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    choices = sample(mut_op.mutation_weights,n_mut,replace = false)

    for index in choices

        proposal = new_w[index] + rand(mut_op.noise_distribution)*new_w[index]

        while (sign(proposal) != sign(new_w[index])) || (abs(proposal) > mut_op.max_w)
            proposal = new_w[index] + rand(mut_op.noise_distribution)*new_w[index]
        end

        new_w[index] = proposal
    end

    valid = maximum(abs.(new_w)) <= mut_op.max_w

    return new_w,nothing,nothing,nothing, valid
end

function noise_add(w::Matrix{Float64},mut_op::MutationOperator)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    choices = sample(mut_op.mutation_weights,n_mut,replace = false)

    for index in choices
        proposal = new_w[index] + rand(mut_op.noise_distribution)
        new_w[index] = abs(proposal) > mut_op.max_w ? mut_op.max_w*sign(proposal) : proposal
    end

    valid = maximum(abs.(new_w)) <= mut_op.max_w

    return new_w,nothing,nothing,nothing, valid
end