using OptimalTransport
using Distances

# Optimal transport approaches

function generate_stripe(positions::Vector{Float64}, widths::Vector{Float64},max_value::Float64)

    tissue = zeros(Nc)

    for (p,w) in zip(positions,widths)
        for i in Int.(p-w/2):Int.(p+w/2)
            tissue[i] = max_value
        end
    end

    tissue
end

function generate_stripe_distribution(positions::Vector{Float64}, widths::Vector{Float64})

    tissue = zeros(Nc)

    for (p,w) in zip(positions,widths)
        for i in Int.(p-w/2):Int.(p+w/2)
            tissue[i] = 1.
        end
    end

    DiscreteNonParametric(1:Nc, tissue ./ sum(tissue))
end


function conc2dist(concentrations::Vector{Float64})
    DiscreteNonParametric(1:length(concentrations), concentrations ./ sum(concentrations))
end

function penalty(x::Vector{Float64},thresh::Float64,n_stripe::Int64,penalisation::Float64)

    up = 0.
    down = 0.

    for i in 1:length(x)-1
        if x[i] <= thresh && x[i+1] > thresh
            up += 1
        elseif x[i] >= thresh && x[i+1] < thresh
            down += 1
        end
    end

    penalisation*(2*n_stripe - up - down)^2 + 1.

end

f_ot_max(target) = x->(1/maximum(x))*ot_cost(SqEuclidean(),conc2dist(x), target)

f_ot_pen(target,thresh,n_stripe,penalisation) = x->penalty(x,thresh,n_stripe,penalisation)*ot_cost(SqEuclidean(),conc2dist(x), target)


# Dual evaluation approach

function f_sim(x::Vector{Float64},thresh::Float64,n_stripe::Int64)

    up = 0.
    down = 0.

    sl = 0.

    segment_lengths = []

    for i in 1:length(x)-1
        sl +=1
        if x[i] <= thresh && x[i+1] > thresh
            up += 1.
            push!(segment_lengths,sl)
            sl = 0.
        elseif x[i] >= thresh && x[i+1] < thresh
            down += 1.
            push!(segment_lengths,sl)
            sl = 0.
        end
    end

    push!(segment_lengths,sl)

    return (-(2*n_stripe - up - down)^2, -sum((segment_lengths .- Nc/(2*n_stripe + 1)) .^ 2)), segment_lengths

end

function f_sim(x::Vector{Float64},thresh::Float64,n_stripe::Int64, target_segment_lengths::Vector{Float64})

    up = 0.
    down = 0.

    sl = 0.

    segment_lengths = []

    for i in 1:length(x)-1
        sl +=1
        if  x[i] <= thresh && x[i+1] >  thresh
            up += 1.
            push!(segment_lengths,sl)
            sl = 0.
        elseif x[i] >= thresh && x[i+1] <  thresh
            down += 1.
            push!(segment_lengths,sl)
            sl = 0.
        end
    end

    push!(segment_lengths,sl)

    if length(segment_lengths) != length(target_segment_lengths)

        return ((2*n_stripe - up - down)^2, sum((segment_lengths .- mean(target_segment_lengths)) .^ 2)), segment_lengths

    else
        return ((2*n_stripe - up - down)^2, sum((segment_lengths .- target_segment_lengths) .^ 2)), segment_lengths
    end

end

function f_sim(x::Vector{Float64},thresh::Float64,n_stripe::Int64, target_segment_lengths::Vector{Float64},min_width::Float64)

    up = 0.
    down = 0.

    sl = 0.

    segment_lengths = []

    for i in 1:length(x)-1
        sl +=1
        if  x[i] <= thresh && x[i+1] > thresh
            up += 1.
            if sl >= min_width
                push!(segment_lengths,sl)
            else
                push!(segment_lengths,0.)
            end
            sl = 0.

        elseif x[i] >= thresh && x[i+1] < thresh
            down += 1.
            if sl >= min_width
                push!(segment_lengths,sl)
            else
                push!(segment_lengths,0.)
            end
            sl = 0.

        end
    end

    push!(segment_lengths,sl)

    if length(segment_lengths) != length(target_segment_lengths)

        return (-(2*n_stripe - up - down)^2, -sum((segment_lengths .- mean(target_segment_lengths)) .^ 2)), segment_lengths

    else
        return (-(2*n_stripe - up - down)^2, -sum((segment_lengths .- target_segment_lengths) .^ 2)), segment_lengths
    end

end

function f_sim_cw(x::Vector{Float64},thresh::Float64,n_stripe::Int64, target_centre_widths::Vector{Tuple{Float64,Float64}},min_width::Float64)

    up = false

    sl = 0.

    high_segment_lengths = []
    low_segment_lengths = []

    for i in 1:length(x)-1
        sl +=1
        if  x[i] <= thresh && x[i+1] > thresh
            up = true
            if sl >= min_width
                push!(low_segment_lengths,sl)
            else
                push!(low_segment_lengths,0.)
            end

            sl = 0.

        elseif x[i] >= thresh && x[i+1] < thresh
            up = false
            if sl >= min_width
                push!(high_segment_lengths,sl)
            else
                push!(high_segment_lengths,0.)
            end

            sl = 0.
        end

    end

    if up
        push!(high_segment_lengths,sl)
    else
        push!(low_segment_lengths,sl)
    end

    n_stripe_pheno = length(high_segment_lengths)
    valid_pheno = length(low_segment_lengths) > 1 ? length(low_segment_lengths) == n_stripe_pheno + 1 : false

    if valid_pheno

        pheno_centre_widths = []

        length_traversed = 0.

        for il in 1:length(low_segment_lengths) - 1

            length_traversed += low_segment_lengths[il]

            found_centre = length_traversed + 0.5*high_segment_lengths[il]

            push!(pheno_centre_widths,(found_centre,0.5*high_segment_lengths[il]))

            length_traversed += high_segment_lengths[il]
        end

        error = 0.

        for (c,w) in target_centre_widths
            error += minimum([(c - cp)^2 + (w - wp)^2 for (cp,wp) in pheno_centre_widths])
        end

        return (-float((n_stripe - n_stripe_pheno)^2), -error), pheno_centre_widths

    else
        # return (-float((n_stripe - n_stripe_pheno)^2), -Inf), [(NaN,NaN)]
        return (-Inf, -Inf), [(NaN,NaN)]
    end

end


function f_sim_eval(x::Vector{Float64},thresh::Float64,min_width ::Float64)

    up = 0.
    down = 0.

    sl = 0.

    segment_lengths = []

    for i in 1:length(x)-1
        sl +=1
        if x[i] <= thresh && x[i+1] > thresh
            up += 1.
            if sl >= min_width
                push!(segment_lengths,sl)
            else
                push!(segment_lengths,0.)
            end
            sl = 0.
        elseif x[i] >= thresh && x[i+1] < thresh
            down += 1.
            if sl >= min_width
                push!(segment_lengths,sl)
            else
                push!(segment_lengths,0.)
            end
            sl = 0.
        end
    end

    push!(segment_lengths,sl)

    return up,down,segment_lengths

end

function f_sim_cw(x::Vector{Float64},thresh::Float64,n_stripe::Int64,target::DiscreteNonParametric, min_width::Float64)

    up = false

    sl = 0.

    high_segment_lengths = []
    low_segment_lengths = []

    for i in 1:length(x)-1
        sl +=1
        if  x[i] <= thresh && x[i+1] > thresh
            up = true
            if sl >= min_width
                push!(low_segment_lengths,sl)
            else
                push!(low_segment_lengths,0.)
            end

            sl = 0.

        elseif x[i] >= thresh && x[i+1] < thresh
            up = false
            if sl >= min_width
                push!(high_segment_lengths,sl)
            else
                push!(high_segment_lengths,0.)
            end

            sl = 0.
        end

    end

    if up
        push!(high_segment_lengths,sl)
    else
        push!(low_segment_lengths,sl)
    end

    n_stripe_pheno = length(high_segment_lengths)
    valid_pheno = length(low_segment_lengths) > 1 ? length(low_segment_lengths) == n_stripe_pheno + 1 : false

    if valid_pheno

        pheno_centre_widths = []

        length_traversed = 0.

        for il in 1:length(low_segment_lengths) - 1

            length_traversed += low_segment_lengths[il]

            found_centre = length_traversed + 0.5*high_segment_lengths[il]

            push!(pheno_centre_widths,(found_centre,0.5*high_segment_lengths[il]))

            length_traversed += high_segment_lengths[il]
        end

        return (-float((n_stripe - n_stripe_pheno)^2), -ot_cost(SqEuclidean(),conc2dist(x), target)), pheno_centre_widths

    else
        # return (-float((n_stripe - n_stripe_pheno)^2), -Inf), [(NaN,NaN)]
        return (-Inf, -Inf), [(NaN,NaN)]
    end

end

function f_mse(x::Vector{Float64},target::Vector{Float64})
    return (0.,-sum((x .- target).^2)), [(NaN,NaN)]

end

function f_dist(x::Vector{Float64},target::DiscreteNonParametric)
    return (0., -ot_cost(SqEuclidean(),conc2dist(x), target)), [(NaN,NaN)]

end


function malt_fitness_half(conc::Vector{Float64})

    Lt = length(conc)

    N = 2

    id_segments = [Int(floor((k-1)*(Lt/N) + 1)) : Int(floor(k*Lt/N)) for k in 1:N]

    high_sum = 0.
    low_sum = 0.

    for i in 1:N
        if i % 2 == 0
            high_sum += sum(conc[id_segments[i]])
        else
            low_sum += sum(conc[id_segments[i]])
        end
    end

    # max_conc = maximum(conc)

    # max_conc = max(min_height,maximum(conc))

    max_conc = 20.

    return (high_sum - low_sum) / ((Lt/N)*max_conc),  [(NaN,NaN)]

end

function malt_fitness(conc::Vector{Float64},n_stripe::Int64)

    Lt = length(conc)

    N = 2*n_stripe + 1

    id_segments = [Int(floor((k-1)*(Lt/N) + 1)) : Int(floor(k*Lt/N)) for k in 1:N]

    high_sum = 0.
    low_sum = 0.

    for i in 1:N
        if i % 2 == 0
            high_sum += sum(conc[id_segments[i]])
        else
            low_sum += sum(conc[id_segments[i]])
        end
    end

    # max_conc = maximum(conc)

    # max_conc = max(min_height,maximum(conc))

    max_conc = 20.

    return ((2/(N-1))*high_sum - (2/(N+1))*low_sum) / ((Lt/N)*max_conc), [(NaN,NaN)]
    # return ((2/(N-1))*high_sum - low_sum) / ((Lt/N)*max_conc), [(NaN,NaN)]

end


function malt_fitness(conc::Vector{Float64},n_stripe::Int64,min_height::Float64)

    Lt = length(conc)

    N = 2*n_stripe + 1

    id_segments = [Int(floor((k-1)*(Lt/N) + 1)) : Int(floor(k*Lt/N)) for k in 1:N]

    high_sum = 0.
    low_sum = 0.

    for i in 1:N
        if i % 2 == 0
            high_sum += sum(conc[id_segments[i]])
        else
            low_sum += sum(conc[id_segments[i]])
        end
    end

    # max_conc = maximum(conc)

    max_conc = max(min_height,maximum(conc))

    return ((2/(N-1))*high_sum - (2/(N+1))*low_sum) / ((Lt/N)*max_conc), [(NaN,NaN)]
    # return ((2/(N-1))*high_sum - low_sum) / ((Lt/N)*max_conc), [(NaN,NaN)]

end

function malt_fitness_right(conc::Vector{Float64})

    Lt = length(conc)

    N = 2

    id_segments = [Int(floor((k-1)*(Lt/N) + 1)) : Int(floor(k*Lt/N)) for k in 1:N]

    high_sum = sum(conc[id_segments[2]])

    low_sum = sum(conc[id_segments[1]])

    # max_conc = maximum(conc)

    # max_conc = max(min_height,maximum(conc))

    max_conc = 20.

    return (high_sum - low_sum) / ((Lt/N)*max_conc),  [(NaN,NaN)]

end

function malt_fitness_left(conc::Vector{Float64})

    Lt = length(conc)

    N = 2

    id_segments = [Int(floor((k-1)*(Lt/N) + 1)) : Int(floor(k*Lt/N)) for k in 1:N]

    high_sum = sum(conc[id_segments[1]])

    low_sum = sum(conc[id_segments[2]])

    # max_conc = maximum(conc)

    # max_conc = max(min_height,maximum(conc))

    max_conc = 20.

    return (high_sum - low_sum) / ((Lt/N)*max_conc),  [(NaN,NaN)]

end


function perfect_malt_conc!(conc::Vector{Float64},n_stripe::Int64,max_conc::Float64)

    Lt = length(conc)

    N = 2*n_stripe + 1

    id_segments = [Int(floor((k-1)*(Lt/N) + 1)) : Int(floor(k*Lt/N)) for k in 1:N]

    high_sum = 0.
    low_sum = 0.

    for i in 1:N
        if i % 2 == 0
            conc[id_segments[i]] .= max_conc
        else
            conc[id_segments[i]] .= 0.
        end
    end
end

function perfect_malt_conc!(conc::Vector{Float64},n_stripe::Int64,max_conc::Float64,k)

    Lt = length(conc)

    N = 2*n_stripe + 1

    id_segments = [Int(floor((k-1)*(Lt/N) + 1)) : Int(floor(k*Lt/N)) for k in 1:N]

    high_sum = 0.
    low_sum = 0.

    for i in 1:N
        if i % k == 0
            conc[id_segments[i]] .= max_conc
        else
            conc[id_segments[i]] .= 0.
        end
    end
end
