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