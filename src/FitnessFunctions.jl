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

    return (2*n_stripe - up - down)^2, sum((segment_lengths .- Nc/(2*n_stripe + 1)) .^ 2), segment_lengths

end

function f_sim(x::Vector{Float64},thresh::Float64,n_stripe::Int64, target_segment_lengths::Vector{Float64})

    up = 0.
    down = 0.

    sl = 0.

    segment_lengths = []

    for i in 1:length(x)-1
        sl +=1
        if  x[i] <= thresh && x[i+1] > thresh
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

    if length(segment_lengths) != length(target_segment_lengths)

        return (2*n_stripe - up - down)^2, sum((segment_lengths .- mean(target_segment_lengths)) .^ 2), segment_lengths

    else
        return (2*n_stripe - up - down)^2, sum((segment_lengths .- target_segment_lengths) .^ 2), segment_lengths
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

        return (2*n_stripe - up - down)^2, sum((segment_lengths .- mean(target_segment_lengths)) .^ 2), segment_lengths

    else
        return (2*n_stripe - up - down)^2, sum((segment_lengths .- target_segment_lengths) .^ 2), segment_lengths
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

function f_sim(x::Vector{Float64},thresh::Float64,n_stripe::Int64,target::DiscreteNonParametric)

    up = 0.
    down = 0.

    for i in 1:length(x)-1
        if x[i] <= thresh && x[i+1] > thresh
            up += 1.
        elseif x[i] >= thresh && x[i+1] < thresh
            down += 1.
        end
    end


    return (2*n_stripe - up - down)^2, ot_cost(SqEuclidean(),conc2dist(x), target)

end
