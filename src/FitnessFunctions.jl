function malt_fitness_relative(conc,n_stripe::Int64)

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

    max_conc = maximum(conc)

    return ((2/(N-1))*high_sum - (2/(N+1))*low_sum) / ((Lt/N)*max_conc)
end

function malt_fitness_absolute(conc,n_stripe::Int64,max_conc::Float64)

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

    return ((2/(N-1))*high_sum - (2/(N+1))*low_sum) / ((Lt/N)*max_conc)
end

function stripe_indicator(conc,min_stripe_width,lower_bound,upper_bound)

    low = findall(conc .< lower_bound)
    high = findall(conc .> upper_bound)

    if (length(low) != 0) & (length(high) != 0)

        valid_arrangment = (low[1] < high[1]) & (low[end] > high[end])
        cts_high_region = !(any([(id > high[1]) & (id < high[end]) for id in low]))
        msw = (high[1] - low[1] >= min_stripe_width) & (low[end] - high[end] >= min_stripe_width) & (high[end] - high[1] + 1 >= min_stripe_width)

        if valid_arrangment & cts_high_region & msw
            return 0.
        else
            return -1.
        end
    else
        return -1. 
    end

end

function add_fitness(tuple_f)
    return tuple_f[1] + ((tuple_f[2]+1)/2) + 1
end

