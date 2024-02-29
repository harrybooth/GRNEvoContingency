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


########################


function noise_specified(w::Vector{Float64},mut_id::Int64,mut_size::Float64,mut_type::Symbol)
    new_w = copy(w)

    if mut_type == :new
        new_w[mut_id] = mut_size
    elseif mut_type == :existing
        new_w[mut_id] = new_w[mut_id] + mut_size*new_w[mut_id]
    elseif mut_type == :del
        new_w[mut_id] = 0
    else
        new_w[mut_id] = NaN
    end

    return new_w
end

function noise_specified(w::Vector{Float64},mut_id::Vector{Int64},mut_size::Vector{Any},mut_type::Vector{Any})

    new_w = copy(w)

    for n in 1:length(mut_id)
        if mut_type[n] == :new
            new_w[mut_id[n]] = mut_size[n]
        elseif mut_type[n] == :existing
            new_w[mut_id[n]] = new_w[mut_id[n]] + mut_size[n]*new_w[mut_id[n]]
        elseif mut_type[n] == :del
            new_w[mut_id[n]] = 0
        else
            new_w[mut_id[n]] = NaN
        end
    end

    return new_w
end

function noise_mtype(w::Matrix{Float64},mut_op::MutationOperator)

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
            n = rand(Uniform(-mut_op.max_w,mut_op.max_w))
            new_w[index] = n
            push!(mtype,:new)
            push!(sizes,n)
        else
            if rand() < mut_op.deletion_p
                push!(sizes,-new_w[index])
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

function noise_mtype_v1(w::Matrix{Float64},mut_op::MutationOperatorNew)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    choices = sample(mut_op.mutation_weights,n_mut,replace = false)
    mtype = []
    sizes = []

    for index in choices
        if new_w[index] == 0
            n = rand(Uniform(-mut_op.max_w,mut_op.max_w))
            new_w[index] = n
            push!(mtype,:new)
            push!(sizes,n)
        else
            n = exp(-1*rand(mut_op.noise_distribution))

            if rand() < mut_op.sign_flip_probability
                new_w[index] = -1*new_w[index]*n
                push!(sizes,-n)
            else
                new_w[index] = new_w[index]*n
                push!(sizes,n)
            end

            if abs(new_w[index]) > mut_op.max_w
                new_w[index] = sign(new_w[index])*mut_op.max_w
            end

            push!(mtype,:existing)
        end
    end

    return new_w, choices, mtype, sizes, true
end

function noise_mtype_v1i(w::Matrix{Float64},mut_op::MutationOperatorNew)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    choices = sample(mut_op.mutation_weights,n_mut,replace = false)
    mtype = []
    sizes = []

    for index in choices
        if rand(Exponential(0.1)) > abs(new_w[index])
            n = rand(Uniform(-mut_op.max_w,mut_op.max_w))
            new_w[index] = n
            push!(mtype,:new)
            push!(sizes,n)
        else
            n = exp(-1*rand(mut_op.noise_distribution))

            if rand() < mut_op.sign_flip_probability
                new_w[index] = -1*new_w[index]*n
                push!(sizes,-n)
            else
                new_w[index] = new_w[index]*n
                push!(sizes,n)
            end

            if abs(new_w[index]) > mut_op.max_w
                new_w[index] = sign(new_w[index])*mut_op.max_w
            end

            push!(mtype,:existing)
        end
    end

    return new_w, choices, mtype, sizes, true
end

function noise_mtype_v1_pm(w::Matrix{Float64},mut_op::MutationOperatorNew)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    choices = sample(mut_op.mutation_weights,n_mut,replace = false)
    mtype = []
    sizes = []

    for index in choices
        if new_w[index] == 0
            n = exp(-1*rand(mut_op.noise_distribution))

            if rand() < mut_op.sign_flip_probability
                new_w[index] = -1*mut_op.start_affinity*n
                push!(sizes,-n)
            else
                new_w[index] = mut_op.start_affinity*n
                push!(sizes,n)
            end

            push!(mtype,:existing)
        else
            n = exp(-1*rand(mut_op.noise_distribution))

            if rand() < mut_op.sign_flip_probability
                new_w[index] = -1*new_w[index]*n
                push!(sizes,-n)
            else
                new_w[index] = new_w[index]*n
                push!(sizes,n)
            end

            if abs(new_w[index]) > mut_op.max_w
                new_w[index] = sign(new_w[index])*mut_op.max_w
            end
        end
    end

    return new_w, choices, mtype, sizes, true
end

function noise_mtype_v3(w::Matrix{Float64},mut_op::MutationOperatorNew)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    choices = sort(sample(mut_op.mutation_weights,n_mut,replace = false))
    mtype = []
    sizes = []

    for index in choices

        if rand() < mut_op.pm_prob

            push!(mtype,:existing)

            if new_w[index] == 0
                n = exp(-1*rand(mut_op.noise_distribution))
                if rand() < 0.5
                    new_w[index] = mut_op.start_affinity*n
                    push!(sizes,n)
                else
                    new_w[index] = -1*mut_op.start_affinity*n
                    push!(sizes,-n)
                end
            else
                n = exp(-1*rand(mut_op.noise_distribution))
                new_w[index] = new_w[index]*n
                push!(sizes,n)
            end

        else
            n = rand(Uniform(-mut_op.max_w,mut_op.max_w))
            new_w[index] = new_w[index] + n
            push!(sizes,n)
            push!(mtype,:new)
        end

    end

    valid = maximum(abs.(new_w)) <= mut_op.max_w

    return new_w, choices, mtype, sizes, valid
end

function noise_mtype_v3_restricted(w::Matrix{Float64},mut_op::MutationOperatorNew)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    choices = sort(sample(mut_op.mutation_weights,n_mut,replace = false))
    mtype = []
    sizes = []

    for index in choices

        if rand() < mut_op.pm_prob

            push!(mtype,:existing)

            if new_w[index] == 0
                n = exp(-1*rand(mut_op.noise_distribution))
                if rand() < 0.5
                    new_w[index] = mut_op.start_affinity*n
                    push!(sizes,n)
                else
                    new_w[index] = -1*mut_op.start_affinity*n
                    push!(sizes,-n)
                end
            else
                n = exp(-1*rand(mut_op.noise_distribution))
                new_w[index] = new_w[index]*n
                push!(sizes,n)
            end

        else
            max_w = maximum(abs.(w))
            n = rand(Uniform(-max_w,max_w))
            new_w[index] = new_w[index] + n
            push!(sizes,n)
            push!(mtype,:new)
        end

    end

    valid = maximum(abs.(new_w)) <= mut_op.max_w

    return new_w, choices, mtype, sizes, valid
end

function noise_mtype_v4(w::Matrix{Float64},mut_op::MutationOperatorNew)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    choices = sort(sample(mut_op.mutation_weights,n_mut,replace = false))
    mtype = []
    sizes = []

    for index in choices

        if rand() < mut_op.pm_prob

            push!(mtype,:existing)

            if new_w[index] == 0
                n = exp(-1*rand(mut_op.noise_distribution))
                if rand() < 0.5
                    new_w[index] = mut_op.start_affinity*n
                    push!(sizes,n)
                else
                    new_w[index] = -1*mut_op.start_affinity*n
                    push!(sizes,-n)
                end
            else
                n = exp(-1*rand(mut_op.noise_distribution))
                new_w[index] = new_w[index]*n
                push!(sizes,n)
            end

        else
            n = rand(Uniform(-mut_op.max_w,mut_op.max_w))
            new_w[index] = n
            push!(sizes,n)
            push!(mtype,:new)
        end

    end

    valid = maximum(abs.(new_w)) <= mut_op.max_w

    return new_w, choices, mtype, sizes, valid
end

function noise_mtype_v4_restricted(w::Matrix{Float64},mut_op::MutationOperatorNew)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    choices = sort(sample(mut_op.mutation_weights,n_mut,replace = false))
    mtype = []
    sizes = []

    for index in choices

        if rand() < mut_op.pm_prob

            push!(mtype,:new)

            if new_w[index] == 0
                n = exp(-1*rand(mut_op.noise_distribution))
                if rand() < 0.5
                    new_w[index] = mut_op.start_affinity*n
                    push!(sizes,n)
                else
                    new_w[index] = -1*mut_op.start_affinity*n
                    push!(sizes,-n)
                end
            else
                n = exp(-1*rand(mut_op.noise_distribution))
                new_w[index] = new_w[index]*n
                push!(sizes,n)
            end

        else
            max_w = maximum(abs.(w))
            n = rand(Uniform(-max_w,max_w))
            new_w[index] = n
            push!(sizes,n)
            push!(mtype,:existing)
        end

    end

    valid = maximum(abs.(new_w)) <= mut_op.max_w

    return new_w, choices, mtype, sizes, valid
end

function noise_mtype_mult_add_2(w::Matrix{Float64},mut_op::MutationOperatorDual)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    mut_op

    choices = sort(sample(mut_op.mutation_weights,n_mut,replace = false))
    mtype = []
    sizes = []

    for index in choices

        if rand() < mut_op.pm_prob

            push!(mtype,:existing)

            if new_w[index] == 0
                n = rand(mut_op.mult_noise_distribution)
                if rand() < 0.5
                    new_w[index] = mut_op.start_affinity*n
                    push!(sizes,n)
                else
                    new_w[index] = -1*mut_op.start_affinity*n
                    push!(sizes,-n)
                end
            else
                n = rand(mut_op.mult_noise_distribution)

                if rand() < mut_op.sign_flip_probability
                    new_w[index] = -1*new_w[index]*n
                    push!(sizes,-n)
                else
                    new_w[index] = new_w[index]*n
                    push!(sizes,n)
                end
            end

        else
            push!(mtype,:new)

            n = rand(mut_op.additive_noise_distribution)

            new_w[index] = abs(new_w[index] + n)
            
            if rand() < mut_op.sign_flip_probability
                new_w[index] = -1*abs(new_w[index] + n)
                push!(sizes,-n)
            else
                new_w[index] = abs(new_w[index] + n)
                push!(sizes,n)
            end

        end

        if abs(new_w[index]) > mut_op.max_w
            new_w[index] = sign(new_w[index])*mut_op.max_w
        end

    end

    return new_w, choices, mtype, sizes, true
end

function noise_mtype_mult_add_3(w::Matrix{Float64},mut_op::MutationOperatorUniform)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    mut_op

    choices = sort(sample(mut_op.mutation_weights,n_mut,replace = false))
    mtype = []
    sizes = []

    for index in choices

        if rand() < mut_op.pm_prob

            push!(mtype,:existing)

            if new_w[index] == 0
                n = rand(mut_op.mult_noise_distribution)
                if rand() < 0.5
                    new_w[index] = mut_op.start_affinity*n
                    push!(sizes,n)
                else
                    new_w[index] = -1*mut_op.start_affinity*n
                    push!(sizes,-n)
                end
            else
                n = rand(mut_op.mult_noise_distribution)
                new_w[index] = new_w[index]*n
                push!(sizes,n)
            end

        else
            push!(mtype,:new)

            n = n = rand(mut_op.additive_noise_distribution)

            new_w[index] = n
            push!(sizes,n)
        end

        if abs(new_w[index]) > mut_op.max_w
            new_w[index] = sign(new_w[index])*mut_op.max_w
        end

    end

    return new_w, choices, mtype, sizes, true
end

function noise_mtype_mult_add_4(w::Matrix{Float64},mut_op::MutationOperatorUniform)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    mut_op

    choices = sort(sample(mut_op.mutation_weights,n_mut,replace = false))
    mtype = []
    sizes = []

    for index in choices

        if rand() < mut_op.pm_prob

            push!(mtype,:existing)

            if new_w[index] == 0
                n = rand(mut_op.mult_noise_distribution)
                if rand() < 0.5
                    new_w[index] = mut_op.start_affinity*n
                    push!(sizes,n)
                else
                    new_w[index] = -1*mut_op.start_affinity*n
                    push!(sizes,-n)
                end
            else
                n = rand(mut_op.mult_noise_distribution)

                if rand() < mut_op.sign_flip_probability
                    new_w[index] = -1*new_w[index]*n
                    push!(sizes,-n)
                else
                    new_w[index] = new_w[index]*n
                    push!(sizes,n)
                end
            end

        else
            push!(mtype,:new)

            n = rand(mut_op.additive_noise_distribution)

            if rand() < mut_op.sign_flip_probability
                new_w[index] = -1*n
                push!(sizes,-n)
            else
                new_w[index] = n
                push!(sizes,n)
            end
        end

        if abs(new_w[index]) > mut_op.max_w
            new_w[index] = sign(new_w[index])*mut_op.max_w
        end

    end

    return new_w, choices, mtype, sizes, true
end

###############


function noise_specified(w::Vector{Float64},mut_id::Vector{Int64},mut_size::Vector{Any},mut_type::Vector{Any},mut_op::MutationOperatorUniform)

    new_w = copy(w)

    for n in 1:length(mut_id)
        if mut_type[n] == :new
            new_w[mut_id[n]] = mut_size[n]
        elseif mut_type[n] == :existing
            if new_w[mut_id[n]] == 0
                new_w[mut_id[n]] = mut_op.start_affinity*mut_size[n]
            else
                new_w[mut_id[n]] = new_w[mut_id[n]]*mut_size[n]
            end
        else
            new_w[mut_id[n]] = NaN
        end

        if abs(new_w[mut_id[n]]) > mut_op.max_w
            new_w[mut_id[n]] = sign(new_w[mut_id[n]])*mut_op.max_w
        end
    end

    return new_w
end

function noise_specified(w::Vector{Float64},mut_id::Vector{Int64},mut_size::Vector{Any},mut_type::Vector{Tuple{Bool, Symbol}},mut_op::MutationOperatorUniform)

    new_w = copy(w)

    for n in 1:length(mut_id)
        if mut_type[n][2] == :additive
            new_w[mut_id[n]] = mut_size[n]
        elseif mut_type[n][2] == :multiplicative
            if new_w[mut_id[n]] == 0
                new_w[mut_id[n]] = mut_op.start_affinity*mut_size[n]
            else
                new_w[mut_id[n]] = new_w[mut_id[n]]*mut_size[n]
            end
        else
            new_w[mut_id[n]] = NaN
        end

        if abs(new_w[mut_id[n]]) > mut_op.max_w
            new_w[mut_id[n]] = sign(new_w[mut_id[n]])*mut_op.max_w
        end
    end

    return new_w
end

function noise_mtype_v1iii(w::Matrix{Float64},mut_op::MutationOperatorNew)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    choices = sample(mut_op.mutation_weights,n_mut,replace = false)
    mtype = []
    sizes = []

    for index in choices

        if new_w[index] == 0
            n = rand(Uniform(-mut_op.max_w,mut_op.max_w))
            new_w[index] = n
            push!(mtype,:new)
            push!(sizes,n)

        elseif rand(Exponential(0.1)) > abs(new_w[index])
            new_w[index] = 0

            push!(mtype,:del)
        else
            n = exp(-1*rand(mut_op.noise_distribution))

            new_w[index] = new_w[index]*n
            push!(sizes,n)
     
            if abs(new_w[index]) > mut_op.max_w
                new_w[index] = sign(new_w[index])*mut_op.max_w
            end

            push!(mtype,:existing)
        end
    end

    return new_w, choices, mtype, sizes, true
end

function noise_mtype_v1iii_sign(w::Matrix{Float64},mut_op::MutationOperatorNew)

    new_w = copy(w)

    n_mut = 0

    while n_mut == 0
        n_mut = mut_op.n_sample_func()
    end

    choices = sample(mut_op.mutation_weights,n_mut,replace = false)
    mtype = []
    sizes = []

    for index in choices

        if new_w[index] == 0
            n = rand(Uniform(-mut_op.max_w,mut_op.max_w))
            new_w[index] = n
            push!(mtype,:new)
            push!(sizes,n)

        elseif rand(Exponential(0.1)) > abs(new_w[index])
            new_w[index] = 0

            push!(mtype,:del)
        else
            n = exp(-1*rand(mut_op.noise_distribution))

            if rand() < mut_op.sign_flip_probability
                new_w[index] = -1*new_w[index]*n
                push!(sizes,-n)
            else
                new_w[index] = new_w[index]*n
                push!(sizes,n)
            end
     
            if abs(new_w[index]) > mut_op.max_w
                new_w[index] = sign(new_w[index])*mut_op.max_w
            end

            push!(mtype,:existing)
        end
    end

    return new_w, choices, mtype, sizes, true
end

function noise_specified(w::Vector{Float64},mut_id::Vector{Int64},mut_size::Vector{Any},mut_type::Vector{Any},mut_op::MutationOperatorNew)

    new_w = copy(w)

    for n in 1:length(mut_id)
        if mut_type[n] == :new
            new_w[mut_id[n]] = mut_size[n]
        elseif mut_type[n] == :existing
            new_w[mut_id[n]] = new_w[mut_id[n]]*mut_size[n]

            if abs(new_w[mut_id[n]]) > mut_op.max_w
                new_w[mut_id[n]] = sign(new_w[mut_id[n]])*mut_op.max_w
            end

        elseif mut_type[n] == :del
            new_w[mut_id[n]] = 0
        else
            new_w[mut_id[n]] = NaN
        end
    end

    return new_w
end

struct MutationOperatorUniform
    mult_noise_distribution :: Distribution
    additive_noise_distribution :: Distribution
    n_sample_func :: Any
    pm_prob :: Float64
    start_affinity :: Float64
    max_w ::Float64
    mutation_weights :: Vector{CartesianIndex{2}}
    sign_flip_probability :: Float64
end

struct MutationOperatorNew
    noise_distribution :: Distribution
    n_sample_func :: Any
    pm_prob :: Float64
    start_affinity :: Float64
    max_w ::Float64
    mutation_weights :: Vector{CartesianIndex{2}}
    sign_flip_probability :: Float64
end