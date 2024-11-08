function set_weight(entry::Tuple{Int,Int},step::Float64,w::Matrix{Float64})
    new_w = copy(w)
    new_w[entry...] = step
    return new_w
end

function set_weight(w,entry_1::Tuple{Int,Int},entry_2::Tuple{Int,Int},weight_1::Float64,weight_2::Float64)
    new_w = copy(w)

    new_w[entry_1...] = weight_1
    new_w[entry_2...] = weight_2

    return vec(new_w)
end

function create_mutant_get_fitness(founder::Individual,development::DESystemSolver,entry_1::Tuple{Int,Int},entry_2::Tuple{Int,Int},weight_1::Float64,weight_2::Float64,fitness_function)

    new_w = copy(founder.genotype.p[1])

    new_w[entry_1...] = weight_1
    new_w[entry_2...] = weight_2

    mutant = Individual(remake(founder.genotype, p = (new_w,founder.genotype.p[2:end]...)),development)

    mutant_fitness = fitness_function(mutant.phenotype)

    return mutant_fitness, vec(new_w)

end

function assign_mst_label(v,top_vertex_map,predict_id,vertex_to_predict_label)
    if isnan(v[1])
        return 0
    else
        if haskey(top_vertex_map,v)
            if top_vertex_map[v] ∈ predict_id
                return vertex_to_predict_label[top_vertex_map[v]]
            else
                return 5
            end
        else
            return 5
        end
    end
end

function get_max_prob_point(n1,n2,pred_grid,prob_grid,fitness_grid,pred_choice,fitness_tol,founder_fitness)

    fit_v = [fixation_probability_kim(f[1] - founder_fitness[1],f[2] - founder_fitness[2],β[1],β[2]) for f in fitness_grid[n1,n2,:]]

    fitness_valid = findall(x->x>fitness_tol,fit_v)
    pred_valid = first.(Tuple.(findall(x->x==pred_choice,pred_grid[n1,n2,:,:])))

    valid = fitness_valid ∩ pred_valid

    if length(valid) > 0
        prob_id = argmax(prob_grid[n1,n2,:,:][valid])

        return valid[prob_id]
    else
        return NaN
    end

end

function get_max_prob_point(n1,n2,pred_grid,prob_grid,fitness_grid,pred_choice,fitness_tol,top_choice,sample_grid_v)

    fit_v = [fixation_probability_kim(f[1] - founder_fitness[1],f[2] - founder_fitness[2],β[1],β[2]) for f in fitness_grid[n1,n2,:]]

    fitness_valid = first.(Tuple.(findall(x->x>fitness_tol,fit_v)))
    pred_valid = first.(Tuple.(findall(x->x==pred_choice,pred_grid[n1,n2,:,:])))
    top_valid_1 = findall(x->sign(x) == sign(top_choice[1]),sample_grid_v[1,:])
    top_valid_2 = findall(x->sign(x) == sign(top_choice[2]),sample_grid_v[2,:])

    valid = fitness_valid ∩ pred_valid ∩ top_valid_1 ∩ top_valid_2

    if length(valid) > 0
        prob_id = argmax(prob_grid[n1,n2,:,:][valid])

        return valid[prob_id]
    else
        print(n1)
        print(n2)
        return NaN
    end

end

function create_pairwise_fitness_landscape(founder,development,model_gtl,test_indices,N_sample,sample_points)

    sample_grid = [[w1,w2] for w1 in sample_points, w2 in sample_points]

    sample_grid_v = reduce(hcat,vec(sample_grid));

    fitness_grid = fill((0.,0.),(10,10,N_sample*N_sample))
    pred_grid  = fill(0,(10,10,N_sample*N_sample))
    entropy_grid  = fill(0.,(10,10,N_sample*N_sample))
    prob_grid = fill(0.,(10,10,N_sample*N_sample))

    for (n1,t1) in enumerate(test_indices)
        for (n2,t2) in enumerate(test_indices)
            if n1 < n2
                zs = [create_mutant_get_fitness(founder,development,t1,t2,w1,w2,fitness_function) for w1 in sample_points, w2 in sample_points];

                netzs = last.(zs)
                fitzs = first.(zs)

                netzs_v = reduce(hcat,vec(netzs))[1:10,:]
                fitzs_v = vec(fitzs)

                gt_dtrain = xgboost.DMatrix(netzs_v |> transpose |> collect, feature_names = weight_names)
                gt_label_probabilities = model_gtl.predict(gt_dtrain)
                gt_label_predictions = mapslices(p->argmax(p),gt_label_probabilities,dims = 2)
                gt_label_entropies = mapslices(p->entropy(p),gt_label_probabilities,dims = 2);

                fitness_grid[n1,n2,:,:] = fitzs_v
                pred_grid[n1,n2,:,:] = gt_label_predictions
                entropy_grid[n1,n2,:,:] = gt_label_entropies
                prob_grid[n1,n2,:,:] = mapslices(p->maximum(p),gt_label_probabilities,dims = 2)
            end
        end
    end

    return fitness_grid,pred_grid,entropy_grid,prob_grid
end

function plot_pairwise_fitness_landscape!(fig,test_indices,sample_points,founder,mut_q,mutation_op,fitness_grid,pred_grid,prob_grid,predict_colors)

    sample_grid = [[w1,w2] for w1 in sample_points, w2 in sample_points]

    sample_grid_v = reduce(hcat,vec(sample_grid));

    start_network = founder.genotype.p[1]
    founder_fitness = fitness_function(founder.phenotype)

    msf = 7.

    all_min_prob = []

    for (n1,t1) in enumerate(test_indices)
        for (n2,t2) in enumerate(test_indices)
            if n1 < n2
                fit_v = log10.([fixation_probability_kim(f[1] - founder_fitness[1],f[2] - founder_fitness[2],β[1],β[2]) for f in fitness_grid[n1,n2,:]])
                # fit_v = log10.([fixation_probability_kim(f - founder_fitness[2],β[1],β[2]) for f in fitness_grid[n1,n2,:]])

                min_prob = minimum(fit_v[fit_v .> log10(2e-53)])
                push!(all_min_prob,min_prob)
            end
        end
    end

    min_prob = minimum(all_min_prob)

    for (n1,t1) in enumerate(test_indices)
        for (n2,t2) in enumerate(test_indices)
            if n1 < n2
                ax = Axis(fig[n2,n1],xlabel = weight_names_latex_m[t1...],ylabel = weight_names_latex_m[t2...])

                sn1 = start_network[t1...]
                sn2 = start_network[t2...]

                # fit_v = log10.([fixation_probability_kim(f - founder_fitness[2],β[1],β[2]) for f in fitness_grid[n1,n2,:]])

                fit_v = log10.([fixation_probability_kim(f[1] - founder_fitness[1],f[2] - founder_fitness[2],β[1],β[2]) for f in fitness_grid[n1,n2,:]])

                pred_vc = [p >  min_prob ? (predict_colors[c[1]],c[2]) : :black for (c,p) in zip(zip(pred_grid[n1,n2,:],prob_grid[n1,n2,:]),fit_v)]

                CairoMakie.scatter!(ax,sample_grid_v, color = pred_vc,markersize = msf)

                CairoMakie.scatter!(ax,reshape([sn1,sn2],(2,1)), color = :yellow, markersize = 2*msf)

                CairoMakie.xlims!(sn1 - quantile(mutation_op.additive_noise_distribution,mut_q),sn1 + quantile(mutation_op.additive_noise_distribution,mut_q))
                CairoMakie.ylims!(sn2 - quantile(mutation_op.additive_noise_distribution,mut_q),sn2 + quantile(mutation_op.additive_noise_distribution,mut_q))

                hidexdecorations!(ax,label = false)
                hideydecorations!(ax,label = false)

                vlines!(ax,0., linestyle = :dash, color = :black)
                hlines!(ax,0., linestyle = :dash, color = :black)

                ax = Axis(fig[n1,n2],xlabel = weight_names_latex_m[t1...],ylabel = weight_names_latex_m[t2...])

                # prob_f = CairoMakie.scatter!(ax,sample_grid_v, color = prob_fv, markersize = msf,colormap = cgrad(:viridis), colorrange = (min_prob,0.), lowclip = :black)
                
                hm = CairoMakie.heatmap!(ax,sample_points,sample_points,reshape(fit_v,(N_sample,N_sample)),colormap = cgrad(:viridis), colorrange = (min_prob,0.), lowclip = :black)
                CairoMakie.scatter!(ax,reshape([sn1,sn2],(2,1)), color = :cyan, markersize = 2*msf)

                if n2 == length(test_indices)
                    cb = Colorbar(fig[:, n2+1], hm;
                    minorticksvisible=true,tickformat = custom_log_formatter,ticklabelsize = 18.
                        )
                end

                CairoMakie.xlims!(ax,sn1 - quantile(mutation_op.additive_noise_distribution,mut_q),sn1 + quantile(mutation_op.additive_noise_distribution,mut_q))
                CairoMakie.ylims!(ax,sn2 - quantile(mutation_op.additive_noise_distribution,mut_q),sn2 + quantile(mutation_op.additive_noise_distribution,mut_q))

                vlines!(ax,sn1, linestyle = :dash, color = :white)
                hlines!(ax,sn2, linestyle = :dash, color = :white)

                hidexdecorations!(ax,label = false)
                hideydecorations!(ax,label = false)
            end
        end
    end

    rowgap!(fig.layout, Relative(0.01))
    colgap!(fig.layout, Relative(0.01))
end

function create_pairwise_mst_fitness_landscape(founder_S0,development,test_indices,N_sample,sample_points)

    s0_fitness = fitness_function(founder_S0.phenotype);

    fitness_grid_S0 = fill((0.,0.),(10,10,N_sample*N_sample))
    stripe_ind_S0 = fill(0,(10,10,N_sample*N_sample))
    mst_grid_S0 = fill(0.,(10,10,N_sample*N_sample,12))

    for (n1,t1) in enumerate(test_indices)
        for (n2,t2) in enumerate(test_indices)
            if n1 < n2

                zs = [create_mutant_get_fitness(founder_S0,development,t1,t2,w1,w2,fitness_function) for w1 in sample_points, w2 in sample_points];

                netzs = last.(zs)
                fitzs = first.(zs)

                fitzs_v = vec(fitzs)
                netzs_v = vec(netzs)

                # stripe_id = map(x->x[1]==1.,fitzs_v)
                stripe_id = map(x->x[1]==0,fitzs_v)

                mst_v = [stripe_id[n] ? select_minimal_topologies(find_minimal_network(vec(net)[1:10],grn_parameters,DefaultGRNSolver(),fitness_function)) : fill(NaN,12) for (n,net) in enumerate(netzs_v)];

                # fitness_grid_S0[n1,n2,:] = map(f->fixation_probability_kim(f[1] - s0_fitness[1],f[2] - s0_fitness[2],β[1],β[2]),fitzs_v)
                fitness_grid_S0[n1,n2,:] = fitzs_v
                stripe_ind_S0[n1,n2,:] = stripe_id

                for (n,net) in enumerate(mst_v)
                    mst_grid_S0[n1,n2,n,:] = net
                end

            end
        end
    end

    return fitness_grid_S0,stripe_ind_S0,mst_grid_S0
end

function plot_mst_fitness_landscape!(fig,test_indices,sample_points,founder,mut_q,mutation_op,fitness_grid,stripe_ind,mst_grid,mst_label_grid,predict_colors)

    sample_grid = [[w1,w2] for w1 in sample_points, w2 in sample_points]

    sample_grid_v = reduce(hcat,vec(sample_grid));

    start_network = founder.genotype.p[1]
    founder_fitness = fitness_function(founder.phenotype)

    msf = 5.

    all_min_prob = []

    for (n1,t1) in enumerate(test_indices)
        for (n2,t2) in enumerate(test_indices)
            if n1 < n2
                # fit_v = log10.([fixation_probability_kim(f[1] - founder_fitness[1],f[2] - founder_fitness[2],β[1],β[2]) for f in fitness_grid[n1,n2,:]])
                fit_v = log10.(fitness_grid[n1,n2,:])

                min_prob = minimum(fit_v[fit_v .> log10(2e-53)])
                push!(all_min_prob,min_prob)
            end
        end
    end

    min_prob = minimum(all_min_prob)


    for (n1,t1) in enumerate(test_indices)
        for (n2,t2) in enumerate(test_indices)
            if n1 < n2

                ax = Axis(fig[n2,n1],xlabel = weight_names_latex_m[t1...],ylabel = weight_names_latex_m[t2...])

                sn1 = start_network[t1...]
                sn2 = start_network[t2...]

                CairoMakie.scatter!(ax,sample_grid_v, color = [n==0 ? :black : predict_colors[Int(n)] for n in mst_label_grid[n1,n2,:]])
                CairoMakie.scatter!(ax,reshape([sn1,sn2],(2,1)), color = :yellow, markersize = 2*msf)

                CairoMakie.xlims!(ax,sn1 - quantile(mutation_op.additive_noise_distribution,mut_q),sn1 + quantile(mutation_op.additive_noise_distribution,mut_q))
                CairoMakie.ylims!(ax,sn2 - quantile(mutation_op.additive_noise_distribution,mut_q),sn2 + quantile(mutation_op.additive_noise_distribution,mut_q))

                hidexdecorations!(ax,label = false)
                hideydecorations!(ax,label = false)

                # vlines!(ax,0., linestyle = :dash, color = :white)
                # hlines!(ax,0., linestyle = :dash, color = :white)

                vlines!(ax,sn1, linestyle = :dash, color = :white)
                hlines!(ax,sn2, linestyle = :dash, color = :white)

                ax = Axis(fig[n1,n2],xlabel = weight_names_latex_m[t1...],ylabel = weight_names_latex_m[t2...])

                # fit_v = log10.([fixation_probability_kim(f[1] - founder_fitness[1],f[2] - founder_fitness[2],β[1],β[2]) for f in fitness_grid[n1,n2,:]])

                fit_v = log10.(fitness_grid[n1,n2,:])

                hm = CairoMakie.heatmap!(ax,sample_points,sample_points,reshape(fit_v,(N_sample,N_sample)),colormap = cgrad(:viridis), colorrange = (min_prob,0.), lowclip = :black)
            
                if n2 == length(test_indices)
                    cb = Colorbar(fig[:, n2+1], hm;
                    minorticksvisible=true,tickformat = custom_log_formatter,ticklabelsize = 18.
                        )
                end

                CairoMakie.scatter!(ax,reshape([sn1,sn2],(2,1)), color = :cyan, markersize = 2*msf)

                CairoMakie.xlims!(ax,sn1 - quantile(mutation_op.additive_noise_distribution,mut_q),sn1 + quantile(mutation_op.additive_noise_distribution,mut_q))
                CairoMakie.ylims!(ax,sn2 - quantile(mutation_op.additive_noise_distribution,mut_q),sn2 + quantile(mutation_op.additive_noise_distribution,mut_q))

                vlines!(ax,sn1, linestyle = :dash, color = :white)
                hlines!(ax,sn2, linestyle = :dash, color = :white)

                hidexdecorations!(ax,label = false)
                hideydecorations!(ax,label = false)
            end
                
        end
    end

    rowgap!(fig.layout, Relative(0.01))
    colgap!(fig.layout, Relative(0.01))
end

function plot_target_network_fitness_landscape!(fig,all_target_networks,top_n,founder,mut_q,mutation_op,fitness_grid)

    pos_fitness = cgrad(:buda)
    neg_fitness = cgrad(:linear_kbc_5_95_c73_n256,rev = true)

    msf = 5.

    fitness_prob_all = []

    beta_choices = zip( [1,100,1,0.01,1] .* β[1] , [1,1,100,1,0.01] .* β[2] )

    beta_exp_names = [L"β=1 \text{ , } N = 10,000",L"β=100 \text{ , } N = 10,000",L"β=1 \text{ , } N = 1,000,000",L"β=0.01 \text{ , } N = 10,000",L"β=1 \text{ , } N = 100"]

    sample_grid = [[w1,w2] for w1 in sample_points, w2 in sample_points]

    sample_grid_v = reduce(hcat,vec(sample_grid));

    founder_fitness = fitness_function(founder.phenotype)

    for n in 1:top_n+1

        start_network = all_target_networks[n][1]
        start_network_m = reshape(start_network,(3,4))

        mut_id = Tuple.(findall(x->x!=0,reshape(all_target_networks[n][2],(3,4)) .!= reshape(all_target_networks[n][1],(3,4))))

        t1 = mut_id[1]
        t2 = mut_id[2]

        n1 = findall(x->x == t1,test_indices)[1]
        n2 = findall(x->x == t2,test_indices)[1]

        prob_all = []

        for βc in beta_choices
            prob_fv = [log10(fixation_probability_kim(0.,f,βc[1],βc[2])) for f in fitness_grid[n1,n2,:] .- founder_fitness[2]]
            push!(prob_all,prob_fv)
        end

        push!(fitness_prob_all,prob_all)
    end

    for n in 1:top_n+1

        start_network = all_target_networks[n][1]
        start_network_m = reshape(start_network,(3,4))

        mut_id = Tuple.(findall(x->x!=0,reshape(all_target_networks[n][2],(3,4)) .!= reshape(all_target_networks[n][1],(3,4))))

        t1 = mut_id[1]
        t2 = mut_id[2]

        n1 = findall(x->x == t1,test_indices)[1]
        n2 = findall(x->x == t2,test_indices)[1]

        sn1 = start_network_m[t1...]
        sn2 = start_network_m[t2...]

        xmin = ceil(Int,sn1 - quantile(mutation_op.additive_noise_distribution,mut_q))
        xmax = ceil(Int,sn1 + quantile(mutation_op.additive_noise_distribution,mut_q))

        ymin = ceil(Int,sn2 - quantile(mutation_op.additive_noise_distribution,mut_q))
        ymax = ceil(Int,sn2 + quantile(mutation_op.additive_noise_distribution,mut_q))

        ###############
        
        axm_fit = Axis(fig[1,n],xlabel = weight_names_latex_m[t1...],ylabel = weight_names_latex_m[t2...])

        fitness_delta = fitness_grid[n1,n2,:] .- founder_fitness[2]

        neg_f = CairoMakie.scatter!(axm_fit,sample_grid_v[:,fitness_delta .< 0], color = [log(abs(f)) for f in fitness_delta[fitness_delta .< 0]], markersize = msf,colormap = neg_fitness)
        pos_f = CairoMakie.scatter!(axm_fit,sample_grid_v[:,fitness_delta .>= 0], color = [f for f in fitness_delta[fitness_delta .>= 0]], markersize = msf,colormap = pos_fitness)

        CairoMakie.scatter!(axm_fit,reshape([sn1,sn2],(2,1)), color = :red, markersize = 10.)

        hidexdecorations!(axm_fit)
        hideydecorations!(axm_fit)

        if n==1
            cb = Colorbar(fig[1, top_n+2], neg_f;
                minorticksvisible=true,width = 10.,ticklabelsize = 12.,tickformat = custom_neglog_formatter,
                    )
            cb = Colorbar(fig[1, top_n+3], pos_f;
                minorticksvisible=true,width = 10.,ticklabelsize = 12.,
                    )
        end

        CairoMakie.xlims!(axm_fit,xmin,xmax)
        CairoMakie.ylims!(axm_fit,ymin,ymax)

        #################

        prob_v = reduce(vcat,fitness_prob_all[n])

        min_prob = minimum(prob_v[prob_v .> log10(2e-53)])

        for i in 1:length(beta_choices)

            if n ==3
                axm_fit = Axis(fig[i+1,n],xlabel = weight_names_latex_m[t1...],ylabel = weight_names_latex_m[t2...],title = beta_exp_names[i])
            else
                axm_fit = Axis(fig[i+1,n],xlabel = weight_names_latex_m[t1...],ylabel = weight_names_latex_m[t2...])
            end

            prob_f = CairoMakie.scatter!(axm_fit,sample_grid_v, color = fitness_prob_all[n][i], markersize = msf,colormap = cgrad(:viridis), colorrange = (min_prob,0.), lowclip = :red)

            if n==1
                cb = Colorbar(fig[i+1, top_n+2], prob_f;
                    minorticksvisible=true,width = 10.,ticklabelsize = 12.,lowclip = :red,tickformat = custom_poslog_formatter,
                        )
            end

            CairoMakie.scatter!(axm_fit,reshape([sn1,sn2],(2,1)), color = :blue, markersize = 10.)

            hidexdecorations!(axm_fit)
            hideydecorations!(axm_fit)

            CairoMakie.xlims!(axm_fit,xmin,xmax)
            CairoMakie.ylims!(axm_fit,ymin,ymax)
        end
    end
end