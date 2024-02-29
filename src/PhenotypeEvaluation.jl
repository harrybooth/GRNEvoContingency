function new_phenotype(ind::Individual,new_w::Matrix{Float64},development)
    Individual(remake(ind.genotype, p = (new_w,ind.genotype.p[2:end]...)),development).phenotype.u[end][3,:]
end

function pheno_characterise(conc,sharp_cell_number,lower_bound,upper_bound,min_stripe_width)

    low = findall(conc .< lower_bound)
    high = findall(conc .> upper_bound)

    left_transition = nothing
    left_sharp = nothing

    right_transition = nothing
    right_sharp = nothing

    max_height = maximum(conc)

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

    valid_pattern = ((length(valid_low)) == (length(valid_high)+1)) & all((conc[1:min_stripe_width] .< lower_bound)) & all(conc[end-min_stripe_width:end] .< lower_bound)

    if (length(low) != 0) & (length(high) != 0)

        if low[1] < high[1] # first low cell is before first high cell --> left boundary
            if high[1] >= min_stripe_width
                left_transition = high[1] # true
                if high[1] - low[1] > sharp_cell_number
                    left_sharp = true
                else
                    left_sharp = false
                end
            end

            if high[end] < low[end]
                if (high[end] < Nc - min_stripe_width) & valid_pattern
                    right_transition = high[end] # true
                    if low[end] - high[end] > sharp_cell_number
                        right_sharp = true
                    end
                end
            end
        else # first high cell is before first low cell --> right boundary, could be left transition but would be inverse
            if low[1] >= min_stripe_width
                right_transition = low[1]
                if low[1] - high[1] > sharp_cell_number
                    right_sharp = true
                else
                    right_sharp = false
                end
            end
        end

    end

    return (left_transition,left_sharp,right_transition,right_sharp,max_height)

end

function characterise_mutation(ph_profile_1,ph_profile_2,thresh_p)

    mutant_profile = [:nothing]

    left_boundary_was_present  = ph_profile_1[1] != nothing
    right_boundary_was_present = ph_profile_1[3] != nothing

    left_boundary_now_present = ph_profile_2[1] != nothing
    right_boundary_now_present = ph_profile_2[3] != nothing


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

        move_boundary = ph_profile_1[3] != ph_profile_2[3]

        if move_boundary
            push!(mutant_profile, :mrb)
        end

    elseif (! right_boundary_now_present) & right_boundary_was_present

        push!(mutant_profile, :drb)
    else
        nothing
    end


    if length(mutant_profile) == 1
        if ph_profile_2[5] > thresh_p*ph_profile_1[5]
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
