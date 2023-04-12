using DataInterpolations
using Memoization
using SparseArrays

mutable struct LocalLandscape
    origin :: Individual    
    origin_fitness :: Union{Float64,Nothing}
    sample_points :: Union{StepRangeLen,Nothing}
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

    origin_fitness,origin_pheno_class = fitness_function(origin.phenotype) 
    
    LocalLandscape(origin,origin_fitness,nothing,nothing,nothing,nothing,nothing)

end

function create_mutant_get_pheno(founder::Individual,development::DESystemSolver,entry::Tuple{Int,Int},step::Float64,fitness_function,noise_application)

    mutant = create_mutant(founder,x->increment_weight(entry,step,x,noise_application),development)

    mutant_fitness,mutant_pheno_class = fitness_function(mutant.phenotype)

    return mutant_fitness

end

function increment_weight(entry::Tuple{Int,Int},step::Float64,w::Matrix{Float64},noise_application)
    new_w = copy(w)
    new_w[entry...] = noise_application(new_w[entry...],step)
    return new_w
end

function compute_slices!(LL::LocalLandscape,range_percentile::Float64,N_sample::Int64,development::DESystemSolver,mutation_op::MutationOperator,fitness_function,noise_application)

    start_stop = quantile.(mutation_op.noise_distribution, [1-range_percentile, range_percentile])

    sample_points = range(start_stop[1],start_stop[2],length = N_sample)

    slice_fitnesses = fill(0.,(size(LL.origin.genotype.p[1])...,N_sample))
    
    @sync for i in 1:size(LL.origin.genotype.p[1],1)
        for j in 1:size(LL.origin.genotype.p[1],2) 
            for s in 1:N_sample
                @spawn slice_fitnesses[i,j,s] = create_mutant_get_pheno(LL.origin,development,(i,j),sample_points[s],fitness_function,noise_application)
            end
        end
    end

    LL.sample_points = sample_points
    LL.slice_fitnesses = slice_fitnesses

end

function calculate_fitness_increase_probability(fitness_slice::Vector{Float64},current_fitness::Float64,sample_points::StepRangeLen,mutation_op::MutationOperator,β::Float64)

    mass = 0.

    dx = step(sample_points)

    for i in 1:length(sample_points)
        Δf = fitness_slice[i] - current_fitness 
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


mutable struct InterpolatedLandscape
    origin :: Individual    
    origin_fitness :: Union{Float64,Nothing}
    n_interp :: Union{Int64,Nothing}
    itp_networks :: Union{Array{Float64, 2},Nothing}
    slice_fitnesses :: Union{Array{Float64, 1},Nothing}
    debug ::Any
end

function InterpolatedLandscape(start_network_v::Vector{Float64},end_network_v::Vector{Float64},n_interp::Int64,grn_parameters::GRNParameters,development::DESystemSolver,fitness_function)

    start_network = reshape(start_network_v,(3,4))

    p = (start_network,grn_parameters.degradation)

    genotype = ODEProblem(gene_regulation_1d!,grn_parameters.g0,(0,Inf),p)
    phenotype  = solve(genotype,development.alg;development.kwargs...)

    origin = Individual(genotype,phenotype)

    origin_fitness,origin_pheno_class = fitness_function(origin.phenotype) 

    ######

    start_end = hcat(start_network_v,end_network_v)

    itp_networks = zeros(size(start_end,1),n_interp)

    for wi in 1:size(start_end,1)
        itp_g = DataInterpolations.LinearInterpolation(start_end[wi,:],[1,n_interp]);
        itp_networks[wi,:] = [itp_g(t) for t in 1:1:n_interp]
    end

    ######

    slice_fitnesses = zeros(n_interp)

    for (n,network) in enumerate(eachcol(itp_networks))

        if n == 1
            slice_fitnesses[n] = origin_fitness
        else
            new_w = reshape(network,(3,4))
            
            mutant = Individual(remake(origin.genotype, p = (new_w,origin.genotype.p[2:end]...)),development)

            if mutant.phenotype.retcode == :Terminated
                mutant_fitness,mutant_pheno_class = fitness_function(mutant.phenotype)
            else
                mutant_fitness = -1.
            end

            slice_fitnesses[n] = mutant_fitness
        end
    end

    InterpolatedLandscape(origin,origin_fitness,n_interp,itp_networks,slice_fitnesses,nothing)

end


function instability(network_1,network_2,grn_parameters,development,fitness_function)

    itp_ll = InterpolatedLandscape(network_1,network_2,30,grn_parameters,development,fitness_function);

    return 0.5*(itp_ll.slice_fitnesses[1] + itp_ll.slice_fitnesses[end] + 2) - minimum(itp_ll.slice_fitnesses .+ 1)

end

mutable struct InterpolatedLandscapeLean
    slice_fitnesses :: Union{Array{Float64, 1},Nothing}
end

function InterpolatedLandscapeLean(start_network_v::Vector{Float64},end_network_v::Vector{Float64},n_interp::Int64,grn_parameters::GRNParameters,development::DESystemSolver,fitness_function)

    start_network = reshape(start_network_v,(3,4))

    p = (start_network,grn_parameters.degradation)

    genotype = ODEProblem(gene_regulation_1d!,grn_parameters.g0,(0,Inf),p)
    phenotype  = solve(genotype,development.alg;development.kwargs...)

    origin = Individual(genotype,phenotype)

    origin_fitness,origin_pheno_class = fitness_function(origin.phenotype) 

    ######

    start_end = hcat(start_network_v,end_network_v)

    itp_networks = zeros(size(start_end,1),n_interp)

    for wi in 1:size(start_end,1)
        itp_g = DataInterpolations.LinearInterpolation(start_end[wi,:],[1,n_interp]);
        itp_networks[wi,:] = [itp_g(t) for t in 1:1:n_interp]
    end

    ######

    slice_fitnesses = zeros(n_interp)

    for (n,network) in enumerate(eachcol(itp_networks))

        if n == 1
            slice_fitnesses[n] = origin_fitness
        else
            new_w = reshape(network,(3,4))
            
            mutant = Individual(remake(origin.genotype, p = (new_w,origin.genotype.p[2:end]...)),development)

            if mutant.phenotype.retcode == :Terminated
                mutant_fitness,mutant_pheno_class = fitness_function(mutant.phenotype)
            else
                mutant_fitness = -1.
            end

            slice_fitnesses[n] = mutant_fitness
        end
    end

    InterpolatedLandscapeLean(slice_fitnesses)

end