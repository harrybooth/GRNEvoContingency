function new_phenotype(ind::Individual,new_w::Matrix{Float64},development)
    Individual(remake(ind.genotype, p = (new_w,ind.genotype.p[2:end]...)),development).phenotype.u[end][3,:]
end

function pheno_characterise(conc,lower_bound,upper_bound,min_stripe_width)

    low = findall(conc .< lower_bound)
    high = findall(conc .> upper_bound)

    left_transition = nothing
    right_transition = nothing

    max_height = round(maximum(conc), digits = 6)

    if (length(low) != 0) & (length(high) != 0)

        if low[1] < high[1] # first low cell is before first high cell --> left boundary
            if high[1] - low[1] >= min_stripe_width
                left_transition = high[1] # true
            end

            if low[end] > high[end]
                if (high[end] - high[1] + 1 >= min_stripe_width) & (low[end] - high[end] >= min_stripe_width)
                    cts_high_region = !(any([(id > high[1]) & (id < high[end]) for id in low]))
                    if cts_high_region
                        right_transition = high[end] # true
                    end
                end
            end
        else # first high cell is before first low cell --> right boundary, could be left transition but would be inverse
            if (low[1] - high[1] >= min_stripe_width)
                right_transition = low[1]
            end
        end

    end

    return (left_transition,right_transition,max_height)

end

function characterise_mutation(ph_profile_1,ph_profile_2,thresh_p)

    mutant_profile = [:neutral]

    left_boundary_was_present  = ph_profile_1[1] != nothing
    right_boundary_was_present = ph_profile_1[2] != nothing

    left_boundary_now_present = ph_profile_2[1] != nothing
    right_boundary_now_present = ph_profile_2[2] != nothing


    if left_boundary_now_present & (! left_boundary_was_present)

        if right_boundary_was_present 
            push!(mutant_profile, :clb_wrb)
        else
            push!(mutant_profile, :clb)
        end

    elseif left_boundary_now_present & left_boundary_was_present

        move_boundary = ph_profile_1[1] != ph_profile_2[1]

        if move_boundary
            push!(mutant_profile, :mlb)
        end

    elseif (! left_boundary_now_present) & left_boundary_was_present

        push!(mutant_profile, :dlb)
    else
        nothing
    end

    #############

    if right_boundary_now_present & (! right_boundary_was_present)

        if left_boundary_was_present 
            push!(mutant_profile, :crb_wlb)
        else
            push!(mutant_profile, :crb)
        end

    elseif right_boundary_now_present & right_boundary_was_present

        move_boundary = ph_profile_1[2] != ph_profile_2[2]

        if move_boundary
            push!(mutant_profile, :mrb)
        end

    elseif (! right_boundary_now_present) & right_boundary_was_present

        push!(mutant_profile, :drb)
    else
        nothing
    end


    if length(mutant_profile) == 1
        if ph_profile_2[3] > thresh_p*ph_profile_1[3]
            push!(mutant_profile,:grow)
        end
    end


    if (:mlb ∈ mutant_profile) & (:mrb ∈ mutant_profile)
        mutant_profile = [:nothing,:mbb]
    end

    if (:clb ∈ mutant_profile) & (:crb ∈ mutant_profile)
        mutant_profile = [:nothing,:cbb]
    end

    if (:dlb ∈ mutant_profile) & (:crb_wlb ∈ mutant_profile)
        mutant_profile = [:nothing,:switch_lr]
    end

    if (:drb ∈ mutant_profile) & (:clb_wrb ∈ mutant_profile)
        mutant_profile = [:nothing,:switch_rl]
    end

    if length(mutant_profile) > 1
        return mutant_profile[2:end]
    else
        return [mutant_profile[1]]
    end
end

