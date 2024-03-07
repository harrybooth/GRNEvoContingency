function reduce_shap_values(tr)

    ftime = tr.weight_edits[1:tr.H0-1] ./ tr.weight_edits[tr.H0]

    n_mut = size(tr.other,1)

    all_imp = []
    all_color = []
    all_weight = []

    if n_mut > 2
        for i in 2:n_mut
            imp = mean(abs.(tr.other[i,:,1:10]) ./ sum(abs.(tr.other[i,:,1:10]),dims = 2),dims = 1)[1,:]
            color = [ftime[i] for _ in 1:10]
            weight_id = [j for j in 1:10]
            push!(all_imp,imp)
            push!(all_color,color)
            push!(all_weight,weight_id)
        end
    else
        for i in 1:n_mut
            imp = mean(abs.(tr.other[i,:,1:10]) ./ sum(abs.(tr.other[i,:,1:10]),dims = 2),dims = 1)[1,:]
            color = [ftime[i] for _ in 1:10]
            weight_id = [j for j in 1:10]
            push!(all_imp,imp)
            push!(all_color,color)
            push!(all_weight,weight_id)
        end
    end

    return reduce(vcat,all_weight),reduce(vcat,all_imp),reduce(vcat,all_color)
end

function shap_values_mut(tr,n_mut)

    ftime = tr.weight_edits[1:tr.H0-1] ./ tr.weight_edits[tr.H0]

    imp = mean(abs.(tr.other[n_mut,:,1:10]) ./ sum(abs.(tr.other[n_mut,:,1:10]),dims = 2),dims = 1)[1,:]
    color = [ftime[n_mut] for _ in 1:10]
    weight_id = [j for j in 1:10]

    return weight_id,imp,color
end

function shap_values_mut_rm(tr,n_mut,label)

    imp = tr.other[n_mut,label,1:10]
    wsign = sign.(tr.geno_traj[n_mut][1:10])
    weight_id = [j for j in 1:10]

    return weight_id,imp,wsign
end

function shap_values_mut_rm(tr,label)

    n_mut = size(tr.other,1)

    all_imp = []
    all_sign = []
    all_weight = []

    for i in 1:n_mut
    
        imp = tr.other[i,label,1:10]
        wsign = sign.(tr.geno_traj[i][1:10])
        weight_id = [j for j in 1:10]

        push!(all_imp,imp)
        push!(all_sign,wsign)
        push!(all_weight,weight_id)
    end

    return reduce(vcat,all_weight),reduce(vcat,all_imp),reduce(vcat,all_sign)
end

label = 2

all_imp = []
all_sign = []
all_weight = []

for tr in filter(x->(x.inc_metagraph_vertices[x.H0] == predict_label_to_vertex[label]) & (vertex_to_predict_label_other[x.gt_label_predictions[2]] == label) & (x.H0 > 2),trajectories_p)
    weight,imp,sign = shap_values_mut_rm(tr,2,label)
    push!(all_imp,imp)
    push!(all_sign,sign)
    push!(all_weight,weight)
end

function get_shap_epi(tr,n_mut)

    tr_mut_weight = tr.mutant_info[n_mut-1].weight_id

    tr_mut_shap = mean(abs.(tr.other[n_mut,:,1:10]) ./ sum(abs.(tr.other[n_mut,:,1:10]),dims = 2),dims = 1)[1,tr_mut_weight]

    return tr_mut_shap,tr.epistasis[n_mut-1][2]
end

function get_shap_epi_symb(tr,epi_symbol)

    mutant_id = findall(x->x[1] == epi_symbol,tr.epistasis[1:tr.H0-2])

    if length(mutant_id) > 0

        all_me = []
        all_ms = []

        for i in mutant_id
            tr_mut_weight = tr.mutant_info[i].weight_id
            tr_mut_shap = mean(abs.(tr.other[i+1,:,1:10]) ./ sum(abs.(tr.other[i+1,:,1:10]),dims = 2),dims = 1)[1,tr_mut_weight]

            push!(all_me,tr.epistasis[i][2])
            push!(all_ms,tr_mut_shap)
        end

        return reduce(vcat,all_ms),reduce(vcat,all_me)
    else
        return [],[]
    end
end

function get_shap_epi_symb(tr,epi_symbol)

    mutant_id = findall(x->x[1] == epi_symbol,tr.epistasis[1:tr.H0-2])

    if length(mutant_id) > 0

        all_me = []
        all_ms = []

        for i in mutant_id
            tr_mut_weight = tr.mutant_info[i].weight_id
            tr_mut_shap = mean(abs.(tr.other[i+1,:,1:10]) ./ sum(abs.(tr.other[i+1,:,1:10]),dims = 2),dims = 1)[1,tr_mut_weight]

            push!(all_me,tr.epistasis[i][2])
            push!(all_ms,tr_mut_shap)
        end

        return reduce(vcat,all_ms),reduce(vcat,all_me)
    else
        return [],[]
    end
end

function get_shap_epi_symb(tr,epi_symbol,label)

    mutant_id = findall(x->x[1] == epi_symbol,tr.epistasis[1:tr.H0-2])

    if length(mutant_id) > 0

        all_me = []
        all_ms = []

        for i in mutant_id
            tr_mut_weight = tr.mutant_info[i].weight_id
            tr_mut_shap = tr.other[i+1,label,tr_mut_weight]

            push!(all_me,tr.epistasis[i][2])
            push!(all_ms,tr_mut_shap)
        end

        return reduce(vcat,all_ms),reduce(vcat,all_me)
    else
        return [],[]
    end
end

# all_imp = reduce(vcat,all_imp)
# all_sign = Int.(reduce(vcat,all_sign) .+ 2)
# all_weight = reduce(vcat,all_weight);

# fig = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
#     resolution = (1000,200),fontsize = 18.)

# ax = Axis(fig[1,1],xgridvisible = false,ygridvisible = false,ylabel = L"\text{Distribution of shapley value magnitudes}")

# CairoMakie.boxplot!(ax,all_weight,all_imp, dodge = all_sign,color = all_sign,show_outliers = false)

# ax.xticks = (1:10,weight_names_latex)

# fig


# all_pred_shap = []
# all_fitness_shap = []
# all_labels = []

# for tr in trajectories
#     for label in 1:4
#         pred_shap,fitness_shap = get_shap_epi_symb(tr,:se,label)
#         push!(all_pred_shap,pred_shap)
#         push!(all_fitness_shap,fitness_shap)
#         push!(all_labels,[label for _ in pred_shap])
#     end
# end

# all_pred_shap = Float64.(reduce(vcat,all_pred_shap))
# all_fitness_shap = Float64.(reduce(vcat,all_fitness_shap));
# all_labels = Int.(reduce(vcat,all_labels));

# grp = all_fitness_shap .< 0.1;

# fig = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
#     resolution = (1000,1000),fontsize = 18.)


# ax = Axis(fig[1,1],xgridvisible = false,ygridvisible = false)

# # CairoMakie.scatter!(ax,all_fitness_shap,all_pred_shap,markersize = 5.)

# # CairoMakie.xlims!(ax,-2,2)

# CairoMakie.boxplot!(ax,grp,abs.(all_pred_shap),show_outliers = false)

# fig