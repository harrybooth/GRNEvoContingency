using DataInterpolations
using Memoization
using SparseArrays
using Distributed

mutable struct LocalLandscape
    origin :: Individual    
    origin_fitness :: Union{Float64,Nothing}
    sample_points :: Any
    slice_fitnesses :: Union{Array{Float64, 3},Nothing}
    slice_phenotypes :: Union{Array{Tuple{Float64,Float64}, 3},Nothing}
    transition_prob :: Union{Array{Float64, 2},Nothing}
    debug ::Any
end


function LocalLandscape(start_network::Matrix{Float64},grn_parameters::GRNParameters,development::DESystemSolver)

    p = (start_network,grn_parameters.degradation)

    genotype = ODEProblem(gene_regulation_1d!,grn_parameters.g0,(0,Inf),p)
    phenotype  = solve(genotype,development.alg;development.kwargs...)

    origin = Individual(genotype,phenotype)

    origin_fitness = fitness_function(origin.phenotype) 
    
    LocalLandscape(origin,origin_fitness,nothing,nothing,nothing,nothing,nothing)

end


function increment_weight(entry::Tuple{Int,Int},step::Float64,w::Matrix{Float64},noise_application)
    new_w = copy(w)
    new_w[entry...] = noise_application(new_w[entry...],step)
    return new_w
end

function set_weight(entry::Tuple{Int,Int},step::Float64,w::Matrix{Float64})
    new_w = copy(w)
    new_w[entry...] = step
    return new_w
end

function create_mutant_get_pheno(founder::Individual,development::DESystemSolver,entry::Tuple{Int,Int},step::Float64,fitness_function,noise_application)

    mutant = create_mutant(founder,x->increment_weight(entry,step,x,noise_application),development)

    mutant_fitness= fitness_function(mutant.phenotype)

    return mutant_fitness

end

# function compute_slices!(LL::LocalLandscape,range_percentile::Float64,N_sample::Int64,development::DESystemSolver,mutation_op::MutationOperator,fitness_function,noise_application)

#     start_stop = quantile.(mutation_op.noise_distribution, [1-range_percentile, range_percentile])

#     sample_points = range(start_stop[1],start_stop[2],length = N_sample)

#     slice_fitnesses = fill(0.,(size(LL.origin.genotype.p[1])...,N_sample))
    
#     @sync for i in 1:size(LL.origin.genotype.p[1],1)
#         for j in 1:size(LL.origin.genotype.p[1],2) 
#             for s in 1:N_sample
#                 @spawn slice_fitnesses[i,j,s] = create_mutant_get_pheno(LL.origin,development,(i,j),sample_points[s],fitness_function,noise_application)
#             end
#         end
#     end

#     LL.sample_points = sample_points
#     LL.slice_fitnesses = slice_fitnesses

# end

function compute_slices!(LL::LocalLandscape,range_percentile::Float64,N_sample::Int64,development::DESystemSolver,mutation_op::MutationOperator,fitness_function,noise_application)

    start_stop = quantile.(mutation_op.noise_distribution, [1-range_percentile, range_percentile])

    # sample_points = range(start_stop[1],start_stop[2],length = N_sample)

    sample_points_1 = range(start_stop[1],0,length = Int(ceil(N_sample/2))) |> collect

    sample_points_2 = range(0,start_stop[2],length = Int(floor(N_sample/2))) |> collect

    sample_points = vcat(sample_points_1,sample_points_2)

    slice_fitnesses = fill(0.,(size(LL.origin.genotype.p[1])...,N_sample))
    
    for i in 1:size(LL.origin.genotype.p[1],1)
        for j in 1:size(LL.origin.genotype.p[1],2) 
            slice_fitnesses[i,j,:] = pmap(sp->create_mutant_get_pheno(LL.origin,development,(i,j),sp,fitness_function,noise_application),sample_points)
        end
    end

    LL.sample_points = sample_points
    LL.slice_fitnesses = slice_fitnesses

end

# function compute_slices!(LL::LocalLandscape,N_sample::Int64,development::DESystemSolver,fitness_function)

#     start_stop = quantile.(mutation_op.noise_distribution, [1-range_percentile, range_percentile])

#     # sample_points = range(start_stop[1],start_stop[2],length = N_sample)

#     sample_points_1 = range(start_stop[1],0,length = Int(ceil(N_sample/2))) |> collect

#     sample_points_2 = range(0,start_stop[2],length = Int(floor(N_sample/2))) |> collect

#     sample_points = vcat(sample_points_1,sample_points_2)

#     slice_fitnesses = fill(0.,(size(LL.origin.genotype.p[1])...,N_sample))
    
#     for i in 1:size(LL.origin.genotype.p[1],1)
#         for j in 1:size(LL.origin.genotype.p[1],2) 
#             slice_fitnesses[i,j,:] = pmap(sp->create_mutant_get_pheno(LL.origin,development,(i,j),sp,fitness_function,noise_application),sample_points)
#         end
#     end

#     LL.sample_points = sample_points
#     LL.slice_fitnesses = slice_fitnesses

# end


function calculate_fitness_increase_probability(fitness_slice::Vector{Float64},current_fitness::Float64,sample_points,mutation_op::MutationOperator,β::Float64)

    mass = 0.

    dx = mean(sample_points[2:end] .- sample_points[1:end-1])

    for i in 1:length(sample_points)
        Δf = max(fitness_slice[i] - current_fitness,0) 
        mass+= (cdf(mutation_op.noise_distribution,sample_points[i] + dx/2) - cdf(mutation_op.noise_distribution,sample_points[i] - dx/2))*fixation_probability(Δf,β)
    end

    return mass
end


function calculate_transition_probabilities!(LL::LocalLandscape,mutation_op::MutationOperator,β::Float64)

    prob = mapslices(s->calculate_fitness_increase_probability(s,LL.origin_fitness,LL.sample_points,mutation_op,β),LL.slice_fitnesses,dims = 3)[:,:,1]

    LL.transition_prob = prob 
end

function calculate_transition_probabilities(LL::LocalLandscape,mutation_op::MutationOperator,β::Float64)

    prob = mapslices(s->calculate_fitness_increase_probability(s,LL.origin_fitness,LL.sample_points,mutation_op,β),LL.slice_fitnesses,dims = 3)[:,:,1]

    return prob
end
    
@memoize Dict function LocalLandscape(start_network::Matrix{Float64},range_percentile::Float64,N_sample::Int,grn_parameters::GRNParameters,development::DESystemSolver,mutation_op::MutationOperator,β::Float64,fitness_function,noise_application)

    print("Calculating Loss Landscape...")
    
    LL = LocalLandscape(start_network,grn_parameters,development)

    compute_slices!(LL,range_percentile,N_sample,development,mutation_op,fitness_function,noise_application)

    calculate_transition_probabilities!(LL,mutation_op,β)

    LL

end

function LocalLandscape_mass(start_network::Matrix{Float64},range_percentile::Float64,N_sample::Int,grn_parameters::GRNParameters,development::DESystemSolver,mutation_op::MutationOperator,β::Float64,fitness_function,noise_application)
    
    LL = LocalLandscape(start_network,grn_parameters,development)

    compute_slices!(LL,range_percentile,N_sample,development,mutation_op,fitness_function,noise_application)

    calculate_transition_probabilities!(LL,mutation_op,β)

    LL

end


