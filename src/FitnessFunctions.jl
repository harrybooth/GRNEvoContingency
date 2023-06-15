function malt_fitness(conc,n_stripe::Int64)

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

    max_conc = 20.

    return ((2/(N-1))*high_sum - (2/(N+1))*low_sum) / ((Lt/N)*max_conc)

end

function malt_fitness(conc,n_stripe::Int64,min_height::Float64)

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

    max_conc = max(min_height,maximum(conc))

    return ((2/(N-1))*high_sum - (2/(N+1))*low_sum) / ((Lt/N)*max_conc)
end

function malt_fitness_right(conc)

    Lt = length(conc)

    N = 2

    id_segments = [Int(floor((k-1)*(Lt/N) + 1)) : Int(floor(k*Lt/N)) for k in 1:N]

    high_sum = sum(conc[id_segments[2]])

    low_sum = sum(conc[id_segments[1]])

    # max_conc = maximum(conc)

    # max_conc = max(min_height,maximum(conc))

    max_conc = 20.

    return (high_sum - low_sum) / ((Lt/N)*max_conc)

end

function malt_fitness_left(conc)

    Lt = length(conc)

    N = 2

    id_segments = [Int(floor((k-1)*(Lt/N) + 1)) : Int(floor(k*Lt/N)) for k in 1:N]

    high_sum = sum(conc[id_segments[1]])

    low_sum = sum(conc[id_segments[2]])

    max_conc = 20.

    return (high_sum - low_sum) / ((Lt/N)*max_conc)

end

function nstripe_fitness(conc,n_stripe::Int64,min_stripe_width,lower_bound,upper_bound)

    low_segments = []
    high_segments = []
    current_low_width = 0.
    current_upper_width = 0.

    for c in conc
        if c < lower_bound
            push!(high_segments,current_upper_width)
            current_low_width += 1.
            current_upper_width = 0.
        elseif c > upper_bound
            push!(low_segments,current_low_width)
            current_low_width = 0.
            current_upper_width += 1.
        end
    end

    push!(high_segments,current_upper_width)
    push!(low_segments,current_low_width)

    valid_low = filter(x->x>=min_stripe_width,low_segments)
    valid_high = filter(x->x>=min_stripe_width,high_segments)

    valid_pattern = (length(valid_low)) == (length(valid_high)+1)

    if valid_pattern
        n_stripe_found = (length(valid_low) + length(valid_high) - 1)/2
    else
        n_stripe_found = 0
    end

    return Float64(-1*abs(n_stripe - n_stripe_found))
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

