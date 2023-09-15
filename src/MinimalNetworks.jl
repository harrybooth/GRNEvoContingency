# function find_minimal_network(network,grn_parameters,development,fitness_function)

#     size_S = length(network)

#     powerset_S = [[bit == '1' ? true : false for bit in string(i;base = 2,pad = size_S)] for i in 0:2^size_S-1]

#     top = []
#     top_sizes = []

#     ind = Individual(reshape(vcat(network,[0.,0.]),(3,4)),grn_parameters,development)

#     for (n,Q) in enumerate(filter(x->x[10] & (x[3] | x[6]),powerset_S))

#         size_Q = sum(Q)

#         new_network = vcat(mask_bool(network,Q),[0.,0.])

#         mutant = Individual(remake(ind.genotype, p = (reshape(new_network,(3,4)),ind.genotype.p[2:end]...)),development)

#         mutant_fitness = fitness_function(mutant.phenotype)
        
#         if (mutant.phenotype.retcode == ReturnCode.Terminated ) & (abs(mutant_fitness[1]) == 0)
#             push!(top,(size_Q,new_network))
#             push!(top_sizes,size_Q)
#         end
#     end

#     min_top = minimum(top_sizes)

#     return top[findall(top_sizes .== min_top)]
# end

function mask(network,topology)

    new_network = copy(network)

    z0 = findall(x->x == 0,topology)

    new_network[z0] .= 0.

    return new_network

end

function mask_bool(network,mask)

    new_network = copy(network)

    z0 = findall(x->!x,mask)

    new_network[z0] .= 0.

    return new_network

end

function find_minimal_network(network,grn_parameters,development,fitness_function)

    size_S = length(network)

    powerset_S = [[bit == '1' ? true : false for bit in string(i;base = 2,pad = size_S)] for i in 0:2^size_S-1]

    top = []
    top_sizes = []

    ind = Individual(reshape(vcat(network,[0.,0.]),(3,4)),grn_parameters,development)

    for (n,Q) in enumerate(filter(x->x[10] & (x[3] | x[6]),powerset_S))

        size_Q = sum(Q)

        new_network = vcat(mask_bool(network,Q),[0.,0.])

        mutant = Individual(remake(ind.genotype, p = (reshape(new_network,(3,4)),ind.genotype.p[2:end]...)),development)

        mutant_fitness = fitness_function(mutant.phenotype)
        
        if SciMLBase.successful_retcode(mutant.phenotype.retcode) && (abs(mutant_fitness[1]) == 0)
            push!(top,(size_Q,new_network))
            push!(top_sizes,size_Q)
        end
    end

    if length(top_sizes) > 0

        min_top = minimum(top_sizes)

        return top[findall(top_sizes .== min_top)]

    else
        return top
    end
end

function mask_by_id(network,keep_id)

    new_network = copy(network)

    z0 = findall(x->!(x ∈ keep_id),1:length(network))

    new_network[z0] .= 0.

    return new_network

end

function evaluate_smj_shapley(smj,start_network,grn_parameters,development,fitness_function)

    non_zero_smj = findall(x->x!=0,smj)

    size_S = length(non_zero_smj)

    powerset_S = [[bit == '1' ? true : false for bit in string(i;base = 2,pad = size_S)] for i in 0:2^size_S-1]

    shapley_values = zeros(length(start_network))

    for (i,w_id) in enumerate(non_zero_smj)

        powerset_S_minus_i = filter(x->!x[i],powerset_S)

        ϕi = 0

        for Q in powerset_S_minus_i 

            size_Q = sum(Q)

            normalize_coefficient = (factorial(size_Q)*factorial(size_S - size_Q - 1))/factorial(size_S)

            Q_add_i = copy(Q)

            Q_add_i[i] = true

            new_network_i = start_network .+ mask_by_id(smj,non_zero_smj[Q_add_i])

            mutant_i = Individual(reshape(new_network_i,(3,4)),grn_parameters,development)

            mutant_fitness_i = fitness_function(mutant_i.phenotype)

            f_Q_add_i = add_fitness(mutant_fitness_i)

            new_network = start_network .+ mask_by_id(smj,non_zero_smj[Q])

            mutant = Individual(reshape(new_network,(3,4)),grn_parameters,development)

            mutant_fitness = fitness_function(mutant.phenotype)

            f_Q = add_fitness(mutant_fitness)
            
            ϕi += normalize_coefficient*(f_Q_add_i - f_Q)
        end

        shapley_values[w_id] = ϕi

    end

    return shapley_values
end