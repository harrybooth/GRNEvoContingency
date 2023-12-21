# https://juliadatascience.io/makie_layouts

function plot_convergence_rate!(ax,conv_time,n_trials,max_gen)

    cum_conv = [sum(conv_time .< i)/n_trials for i in 1:max_gen];

    CairoMakie.lines!(ax,cum_conv,color = :blue,linewidth = 4.)

    ax.title = "A. Convergence Rate"

    ax.xlabel = "Number of generations"
    ax.ylabel = "% of trajectories converged"

end

function plot_mst_distribution!(ax,sorted_counts_uep,top_n,min_n)

    color_sorted_counts_uep = [i <= top_n ? color_scheme[i] : :grey for i in 1:length(sorted_counts_uep)]

    view_sorted_uep_id = sorted_counts_uep .> min_n

    other = sum(sorted_counts_uep[.!view_sorted_uep_id])

    n_other = sum(.!view_sorted_uep_id)

    n_norm = sum(view_sorted_uep_id)

    view_sorted_uep_counts = vcat(sorted_counts_uep[view_sorted_uep_id],[other])

    view_color_sorted_counts_uep = vcat(color_sorted_counts_uep[view_sorted_uep_id],[:grey])

    axis_label = vcat(collect(string.(1:length(view_sorted_uep_counts)-1)),["other"])

    ######################

    CairoMakie.barplot!(ax,view_sorted_uep_counts ./ sum(view_sorted_uep_counts),color = view_color_sorted_counts_uep)

    ax.xticks = (1:length(view_sorted_uep_counts),"MSS-" .* axis_label)

    ax.xticklabelrotation = 45.

    ax.limits = (nothing,(0.,1.))

    CairoMakie.hidedecorations!(ax,label = false,ticklabels = false,ticks = false,minorticks = false)

    ax.xlabel = "Minimal Motif Assignment"
    ax.ylabel = "% of trajectories"
end

function plot_minimal_evo_summary!(fig,trajectories,top_n,sorted_uep,uep_position_dict)

    ax1 = Axis(fig[1,1],title = "A. Convergence Rate",xlabel = "Number of generations")
    ax2 = Axis(fig[1,2],title = "B. Total number of fixed mutations")
    ax3 = Axis(fig[1,3],title = "C. Fitness jumps")

    for (n,parent) in enumerate(sorted_uep[1:top_n])

        tr_data = filter(tr->tr.inc_metagraph_vertices[end] == parent,trajectories)

        nconv = length(tr_data)
    
        conv_time = map(tr->tr.wait_times[end],tr_data)

        cum_conv = [sum(conv_time .< i)/nconv for i in 1:max_gen];

        CairoMakie.lines!(ax1,cum_conv,color = color_scheme[n],linewidth = 4.)
        CairoMakie.lines!(ax1,cum_conv,color = color_scheme[n],linewidth = 4.)
    end

    tr_data = filter(tr->tr.inc_metagraph_vertices[end] ∈ sorted_uep[1:top_n],trajectories)

    n_mutants = map(tr->tr.n_accepted_mutants,tr_data)

    parents = map(tr->uep_position_dict[tr.inc_metagraph_vertices[end]],tr_data)

    grp_colour = map(p->color_scheme[p],parents)

    CairoMakie.boxplot!(ax2,parents,n_mutants,color = grp_colour)

    ax2.xticks = (1:top_n,["MST-" * string(n) for n in 1:top_n])

    fitness_deltas = map(tr->map(x->x[:fitness_delta],tr.mutant_info),tr_data)

    fd_parents = reduce(vcat,[parent.*ones(Int,length(fd)) for (parent,fd) in zip(parents,fitness_deltas)])
    fd_colors = map(p->color_scheme[p],fd_parents)
    fd_data = reduce(vcat,fitness_deltas)

    CairoMakie.violin!(ax3,fd_parents,fd_data,color = fd_colors)

    ax3.xticks = (1:top_n,["MST-" * string(n) for n in 1:top_n])

    for ax in [ax1,ax2,ax3]
        CairoMakie.hidedecorations!(ax,label = false,ticklabels = false,ticks = false,minorticks = false)
    end

end

function draw_example_network!(fig,network,show_pheno,draw_config,node_colors,fs)

    ax = Axis(fig[1,2])
    ax_geno = Axis(fig[1,1],backgroundcolor = RGBf(0.98, 0.98, 0.98))

    example_id = 10

    development = DefaultGRNSolver()

    orig_pheno = Individual(reshape(network,(3,4)),grn_parameters,development).phenotype.u[end]

    if "A" ∈ show_pheno
        CairoMakie.lines!(ax,orig_pheno[1,:], color = node_colors[1], linewidth = 4.)
    end

    if "B" ∈ show_pheno
        CairoMakie.lines!(ax,orig_pheno[2,:], color = node_colors[2], linewidth = 4.)
    end

    if "C" ∈ show_pheno
        CairoMakie.lines!(ax,orig_pheno[3,:], color = node_colors[3], linewidth = 4.)
    end

    if "M" ∈ show_pheno
        morphogen_v = [morph(t) for t in tissue]
        CairoMakie.lines!(ax,morphogen_v, color = node_colors[4], linewidth = 4.)
    end
    
    draw_grn!(ax_geno,orig_network,draw_config,node_colors,fs,false,false)

    hidedecorations!(ax)

end

mutable struct dynamical_summary_config

    fontsize
    embed_markersize
    color_scheme
    node_colors
    draw_config
    fitness_linewidth
    fitness_markersize
    pheno_linewidth
    color_fade
    caption_padding
    caption_fontsize

end

# function dynamical_summary_config(fontsize,embed_markersize,color_scheme,node_colors,draw_config,fitness_linewidth,fitness_markersize,pheno_linewidth,color_fade,caption_padding,caption_fontsize)
#     dynamical_summary_config(fontsize,embed_markersize,color_scheme,node_colors,draw_config,fitness_linewidth,fitness_markersize,pheno_linewidth,color_fade,caption_padding,caption_fontsize)
# end


function plot_dynamical_summary!(fig,trajectories,embedding,top_n,minimal_motif_count,sorted_uep,sorted_counts_uep,mst_conf_int,end_parents,vertex_top_map,example_mst,tr_choice,ds_config)

    trajectories_p_d = filter(tr->tr.inc_metagraph_vertices[end] ∈ sorted_uep[1:top_n],trajectories);

    mo_umap = fig[1:4, 1:4] = GridLayout()

    ex1 = fig[1:3, 5:8] = GridLayout()
    rmh0 = fig[4, 5:8] = GridLayout()

    top_n_dict = Dict(v_id=>pos for (pos,v_id) in enumerate(sorted_uep[1:top_n]))

    #### Motif Distribution

    color_sorted_counts_uep = [i <= top_n ? ds_config.color_scheme[i] : :grey for i in 1:length(sorted_counts_uep)]

    view_sorted_uep_id = sorted_counts_uep .> minimal_motif_count

    other = mean(sorted_counts_uep[.!view_sorted_uep_id])

    n_norm = sum(sorted_counts_uep)

    sorted_uep_proportions = vcat(sorted_counts_uep[view_sorted_uep_id],[other]) ./ n_norm 

    view_color_sorted_uep = vcat(color_sorted_counts_uep[view_sorted_uep_id],[:grey])

    conf_int_choices = mst_conf_int[view_sorted_uep_id]

    push!(conf_int_choices,(minimum(sorted_counts_uep[.!view_sorted_uep_id]) / n_norm,maximum(sorted_counts_uep[.!view_sorted_uep_id]) / n_norm))

    ##############

    # ax1 = Axis(mo_umap[1:2,1:top_n],title = L"\text{Top %$top_n }" * string(top_n) * " MST : " * string(sum(sorted_counts_uep[1:top_n])) * " trajectories", xlabel = L"\text{Dynamics: UMAP 1}", ylabel = L"\text{Dynamics: UMAP 2}")

    count_top_n = round(sum(sorted_uep_proportions[1:top_n])*100, digits = 2)

    ax1 = Axis(mo_umap[1:2,1:top_n],title = L"\text{Top %$top_n  M^{(i)}_{N_i} : %$count_top_n % of trajectories}", xlabel = L"\text{Dynamics: UMAP 1}", ylabel = L"\text{Dynamics: UMAP 2}")

    CairoMakie.scatter!(ax1,embedding, color = [haskey(top_n_dict,i) ? (ds_config.color_scheme[top_n_dict[i]],0.5) : (:grey,0.5) for i in end_parents],markersize = ds_config.embed_markersize)

    hidedecorations!(ax1, label = false)

    for i in 1:top_n

        ax_geno = Axis(mo_umap[3,i], backgroundcolor = (ds_config.color_scheme[i],ds_config.color_fade),aspect = DataAspect())

        draw_grn!(ax_geno,vertex_top_map[sorted_uep[i]],ds_config.draw_config,ds_config.node_colors,ds_config.fontsize,false,false)
    end

    #######################

    ax_mo = Axis(mo_umap[4,1:top_n],ylabel  = L"\text{Probabilty}", xlabel = L"M^{(i)}_{N_i}")

    CairoMakie.barplot!(ax_mo,sorted_uep_proportions,color = view_color_sorted_uep)

    CairoMakie.errorbars!(ax_mo,1:length(sorted_uep_proportions),sorted_uep_proportions,sorted_uep_proportions .- first.(conf_int_choices),last.(conf_int_choices) .- sorted_uep_proportions,color = :black,whiskerwidth = ds_config.fitness_markersize/2)

    ax_mo.xticks = (1:length(sorted_uep_proportions),vcat(string.(1:length(sorted_uep_proportions[1:end-1])),["Other"]))

    CairoMakie.hidedecorations!(ax_mo,label = false,ticklabels = false,ticks = false,minorticks = false)

    ###################

    tr_data = filter(tr->tr.inc_metagraph_vertices[tr.H0] ∈ sorted_uep[example_mst],trajectories);
    tr_data_id = findall(tr->tr.inc_metagraph_vertices[tr.H0] ∈ sorted_uep[example_mst],trajectories)

    tr_traj_id = uniqueid(tr_data[tr_choice].topologies)
    tr_traj_id[end] = length(tr_data[tr_choice].topologies)

    tr_traj = tr_data[tr_choice].topologies[tr_traj_id]
    tr_networks = tr_data[tr_choice].geno_traj[tr_traj_id]

    development = DefaultGRNSolver()

    tr_phenotypes = [Individual(reshape(net,(3,4)),grn_parameters,development).phenotype.u[end] for net in tr_networks]

    tr_fd = create_full_fitness_traj(tr_data[tr_choice].fitness_traj_tuple,tr_data[tr_choice].wait_times)

    tr_top_fitness_id = uniqueid(tr_fd)[tr_traj_id]
    
    tr_fd_coarse = map(x->x[1]+1,tr_fd)
    tr_fd_refine = map(x->x[2],tr_fd)

    tr_top_fitness_rf = [(x,tr_fd_refine[x]) for x in tr_top_fitness_id]

    tr_top_stripe_id = [tr_fd_coarse[x] for x in tr_top_fitness_id]

    ###################

    progression_cs = palette(:haline,length(tr_top_fitness_rf)+1)

    ax_fitness = Axis(ex1[1:2,1:length(tr_phenotypes)],xlabel = L"\text{Generation}")

    hideydecorations!(ax_fitness,label = false,ticklabels = false,ticks = false,minorticks = false)

    rline = CairoMakie.lines!(ax_fitness,tr_fd_refine, color = :grey, linewidth = ds_config.fitness_linewidth)
    cline = CairoMakie.lines!(ax_fitness,tr_fd_coarse, linestyle = "--", color = :blue,linewidth = ds_config.fitness_linewidth)

    CairoMakie.scatter!(ax_fitness,tr_top_fitness_rf, color = [progression_cs[i] for i in 1:length(tr_top_fitness_rf)], markersize = ds_config.fitness_markersize, marker = '★')

    h0 = tr_top_fitness_id[minimum(findall(tr_top_stripe_id .== 1))]
    Ni = tr_top_fitness_id[end]

    if h0 != Ni
        v = Int.(floor((h0+Ni)/2))
        ax_fitness.xticks = ([1,h0,v,Ni],[L"1",L"H_0",L"%$v",L"N_i"])
    else
        ax_fitness.xticks = ([1,h0],[L"1",L"H_0 = N_i"])
    end

    ####################

    # ex1.alignmode = Mixed(right = 0)

    ax_pheno_list = []

    for i in 1:length(tr_phenotypes)

        if tr_top_stripe_id[i] == 1
            ax_geno = Axis(ex1[3:4,i], backgroundcolor = (ds_config.color_scheme[example_mst],ds_config.color_fade),aspect = DataAspect())
        else
            ax_geno = Axis(ex1[3:4,i], backgroundcolor = RGBf(0.98, 0.98, 0.98),aspect = DataAspect())
        end

        ax_pheno = Axis(ex1[5,i],alignmode=Mixed(bottom=0))

        for g in 1:3
            CairoMakie.lines!(ax_pheno,tr_phenotypes[i][g,:],linewidth = ds_config.pheno_linewidth, color = ds_config.node_colors[g])
        end

        CairoMakie.scatter!(ax_pheno,[(90,0.8*tr_phenotypes[end][3,50])],color = progression_cs[i], markersize = ds_config.fitness_markersize,marker = '★')

        CairoMakie.hidedecorations!(ax_pheno)

        draw_grn!(ax_geno,tr_traj[i],ds_config.draw_config,ds_config.node_colors,ds_config.fontsize,false,false)

        push!(ax_pheno_list,ax_pheno)
    end

    linkyaxes!(ax_pheno_list...)

    ##########################

    all_prop  = []
    all_dodge  = []
    all_x  = []

    for n in 1:top_n

        pop = filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories_p_d)

        pop_equal = filter(tr->tr.minimal_stripe_subgraphs[tr.H0] == tr.minimal_stripe_subgraphs[end], pop)

        pop_H0_incl_N = filter(tr->Bool(test_inclusion(tr.minimal_stripe_subgraphs[end],tr.minimal_stripe_subgraphs[tr.H0])) & !(tr.minimal_stripe_subgraphs[end] == tr.minimal_stripe_subgraphs[tr.H0]),pop)

        pop_N_incl_H0 = filter(tr->Bool(test_inclusion(tr.minimal_stripe_subgraphs[tr.H0],tr.minimal_stripe_subgraphs[end])) & !(tr.minimal_stripe_subgraphs[end] == tr.minimal_stripe_subgraphs[tr.H0]),pop)

        n_pop = length(pop)

        proportions = [length(pop_equal),length(pop_H0_incl_N),length(pop_N_incl_H0),length(pop) - length(pop_equal) - length(pop_N_incl_H0) - length(pop_H0_incl_N)]

        @assert sum(proportions) == n_pop

        x = [1,2,3,4]

        dodge = [n,n,n,n]

        push!(all_prop,proportions ./ n_pop)
        push!(all_dodge,dodge)
        push!(all_x,x)

    end

    ax_rh0 = Axis(rmh0[1,1],alignmode=Mixed(top=0))

    x = reduce(vcat,all_x)
    dodge = reduce(vcat,all_dodge)
    proportions = reduce(vcat,all_prop)

    CairoMakie.barplot!(ax_rh0,x,proportions,color = [ds_config.color_scheme[n] for n in dodge],dodge = dodge)

    CairoMakie.hidedecorations!(ax_rh0,label = false,ticklabels = false,ticks = false,minorticks = false)

    ax_rh0.xticks = (1:4,[L"M^{(i)}_{H_{0}} = M^{(i)}_{N_i}",L"M^{(i)}_{H_{0}} \subset M^{(i)}_{N_i}",L"M^{(i)}_{N_i} \subset M^{(i)}_{H_{0}}",L"\text{MST change}"])

    for (label, layout) in zip(["A", "B", "C"], [mo_umap, ex1, rmh0])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = ds_config.caption_fontsize,
            font = :bold,
            padding = (0,ds_config.caption_padding, ds_config.caption_padding, 0),
            halign = :right)
    end
    
    colgap = 5
    rowgap = 10

    colgap!(mo_umap,colgap)
    rowgap!(mo_umap, rowgap)

    colgap!(ex1, colgap)
    rowgap!(ex1, rowgap)

    colgap!(rmh0,colgap)
    rowgap!(rmh0, rowgap)

    rowgap!(fig.layout, Relative(0.01))
    colgap!(fig.layout, Relative(0.01))

end

mutable struct evo_summary_config

    fontsize
    wait_markersize
    wait_linewidth
    color_scheme
    node_colors
    draw_config
    color_fade

    pie_radius
    pie_inner_radius
    pie_colors
    pie_strokewidth
end

function create_evo_summary!(fig,trajectories,top_n,mutation_operator::MutationOperator,sorted_uep, vertex_top_map,wait_time_summary,evo_config)

    # all_wait_times = reduce(hcat,[cumulative_wait_time(tr) for tr in trajectories]);

    all_wait_times = reduce(hcat,[average_wait_time(tr) for tr in trajectories]);

    ax_wait_list = []

    ax_wait_list = []
    ax_wait_2_list = []

    wt_l_list = []
    wt_s_list = []

    min_t_u = -mutation_operator.max_w
    max_t_u = mutation_operator.max_w

    for n in 1:top_n

        # plot_geno = fig[n, 1] = GridLayout()
        # plot_wait = fig[n, 2] = GridLayout()
        # plot_mut_hist = fig[n, 3:4] = GridLayout()
        # plot_epi_types = fig[n, 5:6] = GridLayout()

        plot_geno = fig[n, 1] = GridLayout()
        plot_wait = fig[n, 2:4] = GridLayout()
        plot_mut_hist = fig[n, 5:9] = GridLayout()
        plot_epi_types = fig[n, 10:13] = GridLayout()

        if n==1
            ax_geno = Axis(plot_geno[1,1],backgroundcolor = (evo_config.color_scheme[n],evo_config.color_fade),title =L"M^{(i)}_{N_i}",aspect = DataAspect())
        else
            ax_geno = Axis(plot_geno[1,1],backgroundcolor = (evo_config.color_scheme[n],evo_config.color_fade),aspect = DataAspect())
        end

        # top = Int.(reshape(vertex_top_map[sorted_uep[n]],(3,4)))

        # draw_grn_layout!(ax_geno,top,e_width,vertex_size,arrow_size,arrow_shift,sw,fixed_layout,selfedge_size,node_colors,false)

        draw_grn!(ax_geno,vertex_top_map[sorted_uep[n]],evo_config.draw_config,evo_config.node_colors,evo_config.fontsize,false,false)

        if n == 1
            ax_wait = Axis(plot_wait[1,1], title = L"\mathbb{E}[\text{total weight edits}]",yticklabelsize = 0.8*evo_config.fontsize)
        else
            ax_wait = Axis(plot_wait[1,1],yticklabelsize = 0.8*evo_config.fontsize)
        end

        ax_wait_2 = Axis(plot_wait[1,1], yticklabelcolor = :red, yaxisposition = :right,yscale = log10,yticklabelsize = 0.8*evo_config.fontsize)

        hidespines!(ax_wait_2 )
        hideydecorations!(ax_wait_2,label = false,ticklabels = false,ticks = false,minorticks = false)
        hidexdecorations!(ax_wait_2)

        # #############################

        mut_type_prop_all = []
        mut_type_time_labels = []
        mut_type_labels = []

        # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,1,tr.H0-2), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n]) & (tr.H0-2 > 0),trajectories_p))

        mut_type_prop = map(tr->calculate_mut_type_count(get_mut_type(tr,1,tr.H0-2), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n]) & (tr.H0-2 > 0),trajectories))

        mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

        push!(mut_type_prop_all,mut_type_prop_av)
        push!(mut_type_labels, [1,2,3])
        push!(mut_type_time_labels,[1,1,1])

        # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0-1,tr.H0-1), [:existing,:new,:del]),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n] ,trajectories_p))
        mut_type_prop = map(tr->calculate_mut_type_count(get_mut_type(tr,tr.H0-1,tr.H0-1), [:existing,:new,:del]),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n] ,trajectories))

        mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

        push!(mut_type_prop_all,mut_type_prop_av)
        push!(mut_type_labels, [1,2,3])
        push!(mut_type_time_labels,[2,2,2])

        # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0,length(tr.topologies)-1), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n])  & (tr.H0 < length(tr.topologies)),trajectories_p))

        mut_type_prop = map(tr->calculate_mut_type_count(get_mut_type(tr,tr.H0,length(tr.topologies)-1), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n])  & (tr.H0 < length(tr.topologies)),trajectories))

        mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

        push!(mut_type_prop_all,mut_type_prop_av)
        push!(mut_type_labels, [1,2,3])
        push!(mut_type_time_labels,[3,3,3])

        mut_type_prop_all = reduce(vcat,mut_type_prop_all)
        mut_type_time_labels = reduce(vcat,mut_type_time_labels)
        mut_type_labels = reduce(vcat,mut_type_labels); 

        CairoMakie.barplot!(ax_wait,mut_type_time_labels,mut_type_prop_all,stack = mut_type_labels,color = mut_type_labels)

        if n == top_n
            ax_wait.xticks = (1:3,[L"t<H_{0}",L"t=H_{0}",L"t>H_{0}" ])
        else
            hidexdecorations!(ax_wait)
        end
        
        #format y ticks to latex numbers

        CairoMakie.hidexdecorations!(ax_wait,label = false,ticklabels = false,ticks = false,minorticks = false)
        CairoMakie.hideydecorations!(ax_wait,label = false,ticklabels = false,ticks = false,minorticks = false,grid = false)

        push!(ax_wait_list,ax_wait)

        ############################

        sample_id = findall(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)

        if wait_time_summary == :mean

            mean_wait = mean(all_wait_times[:,sample_id],dims = 2)[:,1]

            std_error_wait = std(all_wait_times[:,sample_id],dims = 2)[:,1] ./ sqrt(length(sample_id))

            mean_wait_type_labels = [1,2,3]

            wt_l = CairoMakie.lines!(ax_wait_2,mean_wait_type_labels,mean_wait,color = :red,linewidth = evo_config.wait_linewidth)
            wt_s = CairoMakie.scatter!(ax_wait_2,mean_wait_type_labels,mean_wait,color = :red,markersize = evo_config.wait_markersize)

            CairoMakie.errorbars!(ax_wait_2,1:length(mean_wait),mean_wait,5 * std_error_wait,color = :red,whiskerwidth = evo_config.wait_markersize/2)

        else

            median_wait_time = mapslices(row->quantile(row, [0.5]),all_wait_times[:,sample_id],dims =2)[:,1]
            lq_wait_time = mapslices(row->quantile(row, [0.25]),all_wait_times[:,sample_id],dims =2)[:,1]
            uq_wait_time = mapslices(row->quantile(row, [0.75]),all_wait_times[:,sample_id],dims =2)[:,1]

            median_wait_type_labels = [1,2,3]

            wt_l = CairoMakie.lines!(ax_wait_2,median_wait_type_labels,median_wait_time,color = :red,linewidth = evo_config.wait_linewidth)
            wt_s = CairoMakie.scatter!(ax_wait_2,median_wait_type_labels,median_wait_time,color = :red,markersize = evo_config.wait_markersize)

            CairoMakie.rangebars!(ax_wait_2,1:length(median_wait_time),lq_wait_time,uq_wait_time,color = :red,whiskerwidth = evo_config.wait_markersize/2)
        end

        push!(ax_wait_2_list,ax_wait_2)

        push!(wt_l_list,wt_l)
        push!(wt_s_list,wt_s)

        #############################

        if n == 1
            ax_epi_lH0 = Axis(plot_epi_types[1,1],title = L"t<H_{0}")
        else
            ax_epi_lH0 = Axis(plot_epi_types[1,1])
        end

        epi_counts = reduce(vcat,map(tr->tr.other[1:tr.H0-2],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_lH0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_lH0)

        if n == 1
            ax_epi_H0 = Axis(plot_epi_types[1,2],title = L"t=H_{0}")
        else
            ax_epi_H0 = Axis(plot_epi_types[1,2])
        end

        epi_counts = map(tr->tr.other[tr.H0-1],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_H0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_H0)

        if n == 1
            ax_epi_uH0 = Axis(plot_epi_types[1,3],title = L"t>H_{0}")
        else
            ax_epi_uH0 = Axis(plot_epi_types[1,3])
        end

        epi_counts = reduce(vcat,map(tr->tr.other[tr.H0:end],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_uH0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_uH0)

        ###############################\

        bins = 50

        for type in [:new,:existing]

            if type == :existing
                mut_noise_dist = mutation_operator.noise_distribution;
            else
                mut_noise_dist = Uniform(-mutation_operator.max_w,mutation_operator.max_w);
            end

            mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

            if type == :existing
                if n==1
                    ax1 = Axis(plot_mut_hist[1,1],title = L"t<H_{0}")
                else
                    ax1 = Axis(plot_mut_hist[1,1])
                end
            else
                ax1 = Axis(plot_mut_hist[2,1])
            end

            if type == :existing
                CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[1])
            else
                CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[2])
            end

            min_t = minimum(mut_size)
            max_t = maximum(mut_size)

            norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];

            CairoMakie.lines!(ax1,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
            CairoMakie.vlines!(ax1,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth)

            mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

            if type == :existing
                if n==1
                    ax2 = Axis(plot_mut_hist[1,2],title = L"t=H_{0}")
                else
                    ax2 = Axis(plot_mut_hist[1,2])
                end
            else
                ax2 = Axis(plot_mut_hist[2,2])
            end

            min_t = minimum(mut_size)
            max_t = maximum(mut_size)

            norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];

            if type == :existing
                CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[1])
            else
                CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[2])
            end

            CairoMakie.lines!(ax2,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
            CairoMakie.vlines!(ax2,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth)

            mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,tr.H0,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

            if type == :existing
                if n==1
                    ax3 = Axis(plot_mut_hist[1,3],title = L"t>H_{0}")
                else
                    ax3 = Axis(plot_mut_hist[1,3])
                end
            else
                ax3 = Axis(plot_mut_hist[2,3])
            end

            min_t = minimum(mut_size)
            max_t = maximum(mut_size)

            norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];

            if type == :existing
                CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[1])
            else
                CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[2])
            end

            CairoMakie.lines!(ax3,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
            CairoMakie.vlines!(ax3,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth)

            # linkxaxes!([ax1,ax2,ax3]...)
            # linkyaxes!([ax1,ax2,ax3]...)

            hidedecorations!(ax1)
            hidedecorations!(ax2)
            hidedecorations!(ax3)
        end

        # colgap!(plot_geno,10)
        # rowgap!(plot_geno,10)
    
        # colgap!(plot_wait, 10)
        # rowgap!(plot_wait, 10)

        # colgap!(plot_mut_hist, 10)
        # rowgap!(plot_mut_hist, 2)

        # colgap!(plot_epi_types, 10)
        # rowgap!(plot_epi_types, 10)

        colgap!(plot_geno,Relative(0.01))
        rowgap!(plot_geno,Relative(0.05))
    
        colgap!(plot_wait, Relative(0.01))
        rowgap!(plot_wait, Relative(0.05))

        colgap!(plot_mut_hist, Relative(0.01))
        rowgap!(plot_mut_hist, Relative(0.05))

        colgap!(plot_epi_types, Relative(0.01))
        rowgap!(plot_epi_types, Relative(0.05))
    end

    # labels = [L"\text{RSE}",L"\text{Sign epistasis}",L"\text{No epistasis}",L"\text{Single mutation}"]

    # Legend(fig[top_n+1, 9:13], [PolyElement(color=c) for c in colors], labels, framevisible=false,nbanks = 2,orientation = :horizontal)

    # colors = palette(:viridis, 3)

    # Legend(fig[top_n+1,2:8],
    #     vcat([[wt_s_list[1], wt_l_list[1]]],[PolyElement(color=c) for c in colors]),
    #     vcat([L"\mathbb{E}[\text{time to accept}]"],[L"\text{weight edits : existing}",L"\text{weight edits : new}",L"\text{weight edits : remove}"]),framevisible=false,nbanks = 2,orientation = :horizontal)

    # labels_wait =  vcat([L"\mathbb{E}[\text{time to accept}]"],[L"\text{weight edits : existing}",L"\text{weight edits : new}",L"\text{weight edits : remove}"])

    # labels_wait =  [L"\mathbb{E}[\text{time to accept}]"]

    if wait_time_summary == :mea
        labels_wait =  [L"\mathbb{E}[\text{total generations}]"]
    else
        labels_wait =  [L"\text{total generations - [25%,50%,75%] quantiles}"]
    end

    labels_mut =  [L"\text{weight edits : existing}",L"\text{weight edits : new}"]

    labels_epi  = [L"\text{RSE}",L"\text{Sign epistasis}",L"\text{No epistasis}",L"\text{Single mutation}"]

    labels = vcat(labels_wait,labels_epi)

    symbol_wait = [[wt_s_list[1], wt_l_list[1]]]

    symbol_mut = [PolyElement(color=c) for c in palette(:viridis, 3)[1:2]]

    symbol_epi = [PolyElement(color=c) for c in evo_config.pie_colors]

    symbol_all = vcat(symbol_wait,symbol_epi)

    # Legend(fig[top_n+1, :], symbol_all, labels, framevisible=false,nbanks = 1,orientation = :horizontal,patchsize = (10, 10), rowgap = 10,colgap = 10)

    # Legend(fig[top_n, 4:13, Bottom()], symbol_all, labels, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = 10,colgap = 10)

    legend_row_gap = 2

    Legend(fig[top_n, 2:4, Bottom()], symbol_wait, labels_wait, framevisible=false,nbanks =1,orientation = :horizontal,patchsize = (10, 10),rowgap = legend_row_gap,colgap = 2,padding=(10.0f0, 10.0f0, 0f0, evo_config.fontsize+1.5*legend_row_gap))

    Legend(fig[top_n, 5:9, Bottom()], symbol_mut, labels_mut, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = legend_row_gap,colgap = 2)
    Legend(fig[top_n,  10:13, Bottom()], symbol_epi, labels_epi, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = legend_row_gap,colgap = 2)

    linkyaxes!(ax_wait_list...)
    linkyaxes!(ax_wait_2_list...)

    rowgap!(fig.layout, Relative(0.01))
    colgap!(fig.layout, Relative(0.01))

end

function create_evo_summary!(fig,trajectories,top_n,mutation_operator::Union{MutationOperatorDual,MutationOperatorUniform},sorted_uep, vertex_top_map,wait_time_summary,evo_config)

    # all_wait_times = reduce(hcat,[cumulative_wait_time(tr) for tr in trajectories]);

    all_wait_times = reduce(hcat,[average_wait_time(tr) for tr in trajectories]);

    ax_wait_list = []

    ax_wait_list = []
    ax_wait_2_list = []

    wt_l_list = []
    wt_s_list = []

    min_t_u = -mutation_operator.max_w
    max_t_u = mutation_operator.max_w

    for n in 1:top_n

        # plot_geno = fig[n, 1] = GridLayout()
        # plot_wait = fig[n, 2] = GridLayout()
        # plot_mut_hist = fig[n, 3:4] = GridLayout()
        # plot_epi_types = fig[n, 5:6] = GridLayout()

        plot_geno = fig[n, 1] = GridLayout()
        plot_wait = fig[n, 2:4] = GridLayout()
        plot_mut_hist = fig[n, 5:9] = GridLayout()
        plot_epi_types = fig[n, 10:13] = GridLayout()

        if n==1
            ax_geno = Axis(plot_geno[1,1],backgroundcolor = (evo_config.color_scheme[n],evo_config.color_fade),title =L"M^{(i)}_{N_i}",aspect = DataAspect())
        else
            ax_geno = Axis(plot_geno[1,1],backgroundcolor = (evo_config.color_scheme[n],evo_config.color_fade),aspect = DataAspect())
        end

        # top = Int.(reshape(vertex_top_map[sorted_uep[n]],(3,4)))

        # draw_grn_layout!(ax_geno,top,e_width,vertex_size,arrow_size,arrow_shift,sw,fixed_layout,selfedge_size,node_colors,false)

        draw_grn!(ax_geno,vertex_top_map[sorted_uep[n]],evo_config.draw_config,evo_config.node_colors,evo_config.fontsize,false,false)

        if n == 1
            ax_wait = Axis(plot_wait[1,1], title = L"\mathbb{E}[\text{total weight edits}]",yticklabelsize = 0.8*evo_config.fontsize)
        else
            ax_wait = Axis(plot_wait[1,1],yticklabelsize = 0.8*evo_config.fontsize)
        end

        ax_wait_2 = Axis(plot_wait[1,1], yticklabelcolor = :red, yaxisposition = :right,yscale = log10,yticklabelsize = 0.8*evo_config.fontsize)

        hidespines!(ax_wait_2 )
        hideydecorations!(ax_wait_2,label = false,ticklabels = false,ticks = false,minorticks = false)
        hidexdecorations!(ax_wait_2)

        # #############################

        mut_type_prop_all = []
        mut_type_time_labels = []
        mut_type_labels = []

        # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,1,tr.H0-2), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n]) & (tr.H0-2 > 0),trajectories_p))

        mut_type_prop = map(tr->calculate_mut_type_count(get_mut_type(tr,1,tr.H0-2), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n]) & (tr.H0-2 > 0),trajectories))

        mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

        push!(mut_type_prop_all,mut_type_prop_av)
        push!(mut_type_labels, [1,2,3])
        push!(mut_type_time_labels,[1,1,1])

        # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0-1,tr.H0-1), [:existing,:new,:del]),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n] ,trajectories_p))
        mut_type_prop = map(tr->calculate_mut_type_count(get_mut_type(tr,tr.H0-1,tr.H0-1), [:existing,:new,:del]),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n] ,trajectories))

        mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

        push!(mut_type_prop_all,mut_type_prop_av)
        push!(mut_type_labels, [1,2,3])
        push!(mut_type_time_labels,[2,2,2])

        # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0,length(tr.topologies)-1), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n])  & (tr.H0 < length(tr.topologies)),trajectories_p))

        mut_type_prop = map(tr->calculate_mut_type_count(get_mut_type(tr,tr.H0,length(tr.topologies)-1), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n])  & (tr.H0 < length(tr.topologies)),trajectories))

        mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

        push!(mut_type_prop_all,mut_type_prop_av)
        push!(mut_type_labels, [1,2,3])
        push!(mut_type_time_labels,[3,3,3])

        mut_type_prop_all = reduce(vcat,mut_type_prop_all)
        mut_type_time_labels = reduce(vcat,mut_type_time_labels)
        mut_type_labels = reduce(vcat,mut_type_labels); 

        CairoMakie.barplot!(ax_wait,mut_type_time_labels,mut_type_prop_all,stack = mut_type_labels,color = mut_type_labels)

        if n == top_n
            ax_wait.xticks = (1:3,[L"t<H_{0}",L"t=H_{0}",L"t>H_{0}" ])
        else
            hidexdecorations!(ax_wait)
        end
        
        #format y ticks to latex numbers

        CairoMakie.hidexdecorations!(ax_wait,label = false,ticklabels = false,ticks = false,minorticks = false)
        CairoMakie.hideydecorations!(ax_wait,label = false,ticklabels = false,ticks = false,minorticks = false,grid = false)

        push!(ax_wait_list,ax_wait)

        ############################ noise_distribut

        sample_id = findall(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)

        if wait_time_summary == :mean

            mean_wait = mean(all_wait_times[:,sample_id],dims = 2)[:,1]

            std_error_wait = std(all_wait_times[:,sample_id],dims = 2)[:,1] ./ sqrt(length(sample_id))

            mean_wait_type_labels = [1,2,3]

            wt_l = CairoMakie.lines!(ax_wait_2,mean_wait_type_labels,mean_wait,color = :red,linewidth = evo_config.wait_linewidth)
            wt_s = CairoMakie.scatter!(ax_wait_2,mean_wait_type_labels,mean_wait,color = :red,markersize = evo_config.wait_markersize)

            CairoMakie.errorbars!(ax_wait_2,1:length(mean_wait),mean_wait,5 * std_error_wait,color = :red,whiskerwidth = evo_config.wait_markersize/2)

        else

            median_wait_time = mapslices(row->quantile(row, [0.5]),all_wait_times[:,sample_id],dims =2)[:,1]
            lq_wait_time = mapslices(row->quantile(row, [0.25]),all_wait_times[:,sample_id],dims =2)[:,1]
            uq_wait_time = mapslices(row->quantile(row, [0.75]),all_wait_times[:,sample_id],dims =2)[:,1]

            median_wait_type_labels = [1,2,3]

            wt_l = CairoMakie.lines!(ax_wait_2,median_wait_type_labels,median_wait_time,color = :red,linewidth = evo_config.wait_linewidth)
            wt_s = CairoMakie.scatter!(ax_wait_2,median_wait_type_labels,median_wait_time,color = :red,markersize = evo_config.wait_markersize)

            CairoMakie.rangebars!(ax_wait_2,1:length(median_wait_time),lq_wait_time,uq_wait_time,color = :red,whiskerwidth = evo_config.wait_markersize/2)
        end

        push!(ax_wait_2_list,ax_wait_2)

        push!(wt_l_list,wt_l)
        push!(wt_s_list,wt_s)

        #############################

        if n == 1
            ax_epi_lH0 = Axis(plot_epi_types[1,1],title = L"t<H_{0}")
        else
            ax_epi_lH0 = Axis(plot_epi_types[1,1])
        end

        epi_counts = reduce(vcat,map(tr->tr.other[1:tr.H0-2],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_lH0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_lH0)

        if n == 1
            ax_epi_H0 = Axis(plot_epi_types[1,2],title = L"t=H_{0}")
        else
            ax_epi_H0 = Axis(plot_epi_types[1,2])
        end

        epi_counts = map(tr->tr.other[tr.H0-1],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_H0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_H0)

        if n == 1
            ax_epi_uH0 = Axis(plot_epi_types[1,3],title = L"t>H_{0}")
        else
            ax_epi_uH0 = Axis(plot_epi_types[1,3])
        end

        epi_counts = reduce(vcat,map(tr->tr.other[tr.H0:end],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_uH0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_uH0)

        ###############################\

        bins = 50

        for type in [:new,:existing]

            if type == :existing
                mut_noise_dist = mutation_operator.mult_noise_distribution;
            else
                mut_noise_dist = mutation_operator.additive_noise_distribution;
            end

            mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

            if type == :existing
                if n==1
                    ax1 = Axis(plot_mut_hist[1,1],title = L"t<H_{0}")
                else
                    ax1 = Axis(plot_mut_hist[1,1])
                end
            else
                ax1 = Axis(plot_mut_hist[2,1])
            end

            if type == :existing
                CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[1])
            else
                CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[2])
            end

            min_t = minimum(mut_size)
            max_t = maximum(mut_size)

            norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];

            CairoMakie.lines!(ax1,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
            CairoMakie.vlines!(ax1,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth)

            mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

            if type == :existing
                if n==1
                    ax2 = Axis(plot_mut_hist[1,2],title = L"t=H_{0}")
                else
                    ax2 = Axis(plot_mut_hist[1,2])
                end
            else
                ax2 = Axis(plot_mut_hist[2,2])
            end

            min_t = minimum(mut_size)
            max_t = maximum(mut_size)

            norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];

            if type == :existing
                CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[1])
            else
                CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[2])
            end

            CairoMakie.lines!(ax2,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
            CairoMakie.vlines!(ax2,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth)

            mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,tr.H0,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

            if type == :existing
                if n==1
                    ax3 = Axis(plot_mut_hist[1,3],title = L"t>H_{0}")
                else
                    ax3 = Axis(plot_mut_hist[1,3])
                end
            else
                ax3 = Axis(plot_mut_hist[2,3])
            end

            min_t = minimum(mut_size)
            max_t = maximum(mut_size)

            norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];

            if type == :existing
                CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[1])
            else
                CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[2])
            end

            CairoMakie.lines!(ax3,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
            CairoMakie.vlines!(ax3,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth)

            # linkxaxes!([ax1,ax2,ax3]...)
            # linkyaxes!([ax1,ax2,ax3]...)

            hidedecorations!(ax1)
            hidedecorations!(ax2)
            hidedecorations!(ax3)
        end

        # colgap!(plot_geno,10)
        # rowgap!(plot_geno,10)
    
        # colgap!(plot_wait, 10)
        # rowgap!(plot_wait, 10)

        # colgap!(plot_mut_hist, 10)
        # rowgap!(plot_mut_hist, 2)

        # colgap!(plot_epi_types, 10)
        # rowgap!(plot_epi_types, 10)

        colgap!(plot_geno,Relative(0.01))
        rowgap!(plot_geno,Relative(0.05))
    
        colgap!(plot_wait, Relative(0.01))
        rowgap!(plot_wait, Relative(0.05))

        colgap!(plot_mut_hist, Relative(0.01))
        rowgap!(plot_mut_hist, Relative(0.05))

        colgap!(plot_epi_types, Relative(0.01))
        rowgap!(plot_epi_types, Relative(0.05))
    end

    # labels = [L"\text{RSE}",L"\text{Sign epistasis}",L"\text{No epistasis}",L"\text{Single mutation}"]

    # Legend(fig[top_n+1, 9:13], [PolyElement(color=c) for c in colors], labels, framevisible=false,nbanks = 2,orientation = :horizontal)

    # colors = palette(:viridis, 3)

    # Legend(fig[top_n+1,2:8],
    #     vcat([[wt_s_list[1], wt_l_list[1]]],[PolyElement(color=c) for c in colors]),
    #     vcat([L"\mathbb{E}[\text{time to accept}]"],[L"\text{weight edits : existing}",L"\text{weight edits : new}",L"\text{weight edits : remove}"]),framevisible=false,nbanks = 2,orientation = :horizontal)

    # labels_wait =  vcat([L"\mathbb{E}[\text{time to accept}]"],[L"\text{weight edits : existing}",L"\text{weight edits : new}",L"\text{weight edits : remove}"])

    # labels_wait =  [L"\mathbb{E}[\text{time to accept}]"]

    if wait_time_summary == :mea
        labels_wait =  [L"\mathbb{E}[\text{total generations}]"]
    else
        labels_wait =  [L"\text{total generations - [25%,50%,75%] quantiles}"]
    end

    labels_mut =  [L"\text{weight edits: multiplicative}",L"\text{weight edits: additive}"]

    labels_epi  = [L"\text{RSE}",L"\text{Sign epistasis}",L"\text{No epistasis}",L"\text{Single mutation}"]

    labels = vcat(labels_wait,labels_epi)

    symbol_wait = [[wt_s_list[1], wt_l_list[1]]]

    symbol_mut = [PolyElement(color=c) for c in palette(:viridis, 3)[1:2]]

    symbol_epi = [PolyElement(color=c) for c in evo_config.pie_colors]

    symbol_all = vcat(symbol_wait,symbol_epi)

    # Legend(fig[top_n+1, :], symbol_all, labels, framevisible=false,nbanks = 1,orientation = :horizontal,patchsize = (10, 10), rowgap = 10,colgap = 10)

    # Legend(fig[top_n, 4:13, Bottom()], symbol_all, labels, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = 10,colgap = 10)

    legend_row_gap = 2

    Legend(fig[top_n, 2:4, Bottom()], symbol_wait, labels_wait, framevisible=false,nbanks =1,orientation = :horizontal,patchsize = (10, 10),rowgap = legend_row_gap,colgap = 2,padding=(10.0f0, 10.0f0, 0f0, evo_config.fontsize+1.5*legend_row_gap))

    Legend(fig[top_n, 5:9, Bottom()], symbol_mut, labels_mut, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = legend_row_gap,colgap = 2)
    Legend(fig[top_n,  10:13, Bottom()], symbol_epi, labels_epi, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = legend_row_gap,colgap = 2)

    linkyaxes!(ax_wait_list...)
    linkyaxes!(ax_wait_2_list...)

    rowgap!(fig.layout, Relative(0.01))
    colgap!(fig.layout, Relative(0.01))

end

function create_evo_summary!(fig,trajectories,top_n,mutation_operator::MutationOperatorNew,sorted_uep, vertex_top_map,wait_time_summary,evo_config)

    # all_wait_times = reduce(hcat,[cumulative_wait_time(tr) for tr in trajectories]);

    all_wait_times = reduce(hcat,[average_wait_time(tr) for tr in trajectories]);

    ax_wait_list = []

    ax_wait_list = []
    ax_wait_2_list = []

    wt_l_list = []
    wt_s_list = []

    min_t_u = -mutation_operator.max_w
    max_t_u = mutation_operator.max_w

    for n in 1:top_n

        # plot_geno = fig[n, 1] = GridLayout()
        # plot_wait = fig[n, 2] = GridLayout()
        # plot_mut_hist = fig[n, 3:4] = GridLayout()
        # plot_epi_types = fig[n, 5:6] = GridLayout()

        plot_geno = fig[n, 1] = GridLayout()
        plot_wait = fig[n, 2:4] = GridLayout()
        plot_mut_hist = fig[n, 5:9] = GridLayout()
        plot_epi_types = fig[n, 10:13] = GridLayout()

        if n==1
            ax_geno = Axis(plot_geno[1,1],backgroundcolor = (evo_config.color_scheme[n],evo_config.color_fade),title =L"M^{(i)}_{N_i}",aspect = DataAspect())
        else
            ax_geno = Axis(plot_geno[1,1],backgroundcolor = (evo_config.color_scheme[n],evo_config.color_fade),aspect = DataAspect())
        end

        # top = Int.(reshape(vertex_top_map[sorted_uep[n]],(3,4)))

        # draw_grn_layout!(ax_geno,top,e_width,vertex_size,arrow_size,arrow_shift,sw,fixed_layout,selfedge_size,node_colors,false)

        draw_grn!(ax_geno,vertex_top_map[sorted_uep[n]],evo_config.draw_config,evo_config.node_colors,evo_config.fontsize,false,false)

        if n == 1
            ax_wait = Axis(plot_wait[1,1], title = L"\mathbb{E}[\text{total weight edits}]",yticklabelsize = 0.8*evo_config.fontsize)
        else
            ax_wait = Axis(plot_wait[1,1],yticklabelsize = 0.8*evo_config.fontsize)
        end

        ax_wait_2 = Axis(plot_wait[1,1], yticklabelcolor = :red, yaxisposition = :right,yscale = log10,yticklabelsize = 0.8*evo_config.fontsize)

        hidespines!(ax_wait_2 )
        hideydecorations!(ax_wait_2,label = false,ticklabels = false,ticks = false,minorticks = false)
        hidexdecorations!(ax_wait_2)

        # #############################

        mut_type_prop_all = []
        mut_type_time_labels = []
        mut_type_labels = []

        # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,1,tr.H0-2), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n]) & (tr.H0-2 > 0),trajectories_p))

        mut_type_prop = map(tr->calculate_mut_type_count(get_mut_type(tr,1,tr.H0-2), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n]) & (tr.H0-2 > 0),trajectories))

        mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

        push!(mut_type_prop_all,mut_type_prop_av)
        push!(mut_type_labels, [1,2,3])
        push!(mut_type_time_labels,[1,1,1])

        # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0-1,tr.H0-1), [:existing,:new,:del]),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n] ,trajectories_p))
        mut_type_prop = map(tr->calculate_mut_type_count(get_mut_type(tr,tr.H0-1,tr.H0-1), [:existing,:new,:del]),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n] ,trajectories))

        mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

        push!(mut_type_prop_all,mut_type_prop_av)
        push!(mut_type_labels, [1,2,3])
        push!(mut_type_time_labels,[2,2,2])

        # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0,length(tr.topologies)-1), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n])  & (tr.H0 < length(tr.topologies)),trajectories_p))

        mut_type_prop = map(tr->calculate_mut_type_count(get_mut_type(tr,tr.H0,length(tr.topologies)-1), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n])  & (tr.H0 < length(tr.topologies)),trajectories))

        mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

        push!(mut_type_prop_all,mut_type_prop_av)
        push!(mut_type_labels, [1,2,3])
        push!(mut_type_time_labels,[3,3,3])

        mut_type_prop_all = reduce(vcat,mut_type_prop_all)
        mut_type_time_labels = reduce(vcat,mut_type_time_labels)
        mut_type_labels = reduce(vcat,mut_type_labels); 

        CairoMakie.barplot!(ax_wait,mut_type_time_labels,mut_type_prop_all,stack = mut_type_labels,color = mut_type_labels)

        if n == top_n
            ax_wait.xticks = (1:3,[L"t<H_{0}",L"t=H_{0}",L"t>H_{0}" ])
        else
            hidexdecorations!(ax_wait)
        end
        
        #format y ticks to latex numbers

        CairoMakie.hidexdecorations!(ax_wait,label = false,ticklabels = false,ticks = false,minorticks = false)
        CairoMakie.hideydecorations!(ax_wait,label = false,ticklabels = false,ticks = false,minorticks = false,grid = false)

        push!(ax_wait_list,ax_wait)

        ############################

        sample_id = findall(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)

        if wait_time_summary == :mean

            mean_wait = mean(all_wait_times[:,sample_id],dims = 2)[:,1]

            std_error_wait = std(all_wait_times[:,sample_id],dims = 2)[:,1] ./ sqrt(length(sample_id))

            mean_wait_type_labels = [1,2,3]

            wt_l = CairoMakie.lines!(ax_wait_2,mean_wait_type_labels,mean_wait,color = :red,linewidth = evo_config.wait_linewidth)
            wt_s = CairoMakie.scatter!(ax_wait_2,mean_wait_type_labels,mean_wait,color = :red,markersize = evo_config.wait_markersize)

            CairoMakie.errorbars!(ax_wait_2,1:length(mean_wait),mean_wait,5 * std_error_wait,color = :red,whiskerwidth = evo_config.wait_markersize/2)

        else

            median_wait_time = mapslices(row->quantile(row, [0.5]),all_wait_times[:,sample_id],dims =2)[:,1]
            lq_wait_time = mapslices(row->quantile(row, [0.25]),all_wait_times[:,sample_id],dims =2)[:,1]
            uq_wait_time = mapslices(row->quantile(row, [0.75]),all_wait_times[:,sample_id],dims =2)[:,1]

            median_wait_type_labels = [1,2,3]

            wt_l = CairoMakie.lines!(ax_wait_2,median_wait_type_labels,median_wait_time,color = :red,linewidth = evo_config.wait_linewidth)
            wt_s = CairoMakie.scatter!(ax_wait_2,median_wait_type_labels,median_wait_time,color = :red,markersize = evo_config.wait_markersize)

            CairoMakie.rangebars!(ax_wait_2,1:length(median_wait_time),lq_wait_time,uq_wait_time,color = :red,whiskerwidth = evo_config.wait_markersize/2)
        end

        push!(ax_wait_2_list,ax_wait_2)

        push!(wt_l_list,wt_l)
        push!(wt_s_list,wt_s)

        #############################

        if n == 1
            ax_epi_lH0 = Axis(plot_epi_types[1,1],title = L"t<H_{0}")
        else
            ax_epi_lH0 = Axis(plot_epi_types[1,1])
        end

        epi_counts = reduce(vcat,map(tr->tr.other[1:tr.H0-2],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_lH0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_lH0)

        if n == 1
            ax_epi_H0 = Axis(plot_epi_types[1,2],title = L"t=H_{0}")
        else
            ax_epi_H0 = Axis(plot_epi_types[1,2])
        end

        epi_counts = map(tr->tr.other[tr.H0-1],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_H0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_H0)

        if n == 1
            ax_epi_uH0 = Axis(plot_epi_types[1,3],title = L"t>H_{0}")
        else
            ax_epi_uH0 = Axis(plot_epi_types[1,3])
        end

        epi_counts = reduce(vcat,map(tr->tr.other[tr.H0:end],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_uH0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_uH0)

        ###############################\

        bins = 50

        for type in [:new,:existing]

            if type == :existing
                mut_noise_dist = mutation_operator.noise_distribution;
            else
                mut_noise_dist = Uniform(-mutation_operator.max_w,mutation_operator.max_w);
            end

            mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

            if type == :existing
                if n==1
                    ax1 = Axis(plot_mut_hist[1,1],title = L"t<H_{0}")
                else
                    ax1 = Axis(plot_mut_hist[1,1])
                end
            else
                ax1 = Axis(plot_mut_hist[2,1])
            end

            if type == :existing
                CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[1])
            else
                CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[2])
            end

            min_t = minimum(mut_size)
            max_t = maximum(mut_size)

            norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];

            CairoMakie.lines!(ax1,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
            CairoMakie.vlines!(ax1,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth)

            mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

            if type == :existing
                if n==1
                    ax2 = Axis(plot_mut_hist[1,2],title = L"t=H_{0}")
                else
                    ax2 = Axis(plot_mut_hist[1,2])
                end
            else
                ax2 = Axis(plot_mut_hist[2,2])
            end

            min_t = minimum(mut_size)
            max_t = maximum(mut_size)

            norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];

            if type == :existing
                CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[1])
            else
                CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[2])
            end

            CairoMakie.lines!(ax2,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
            CairoMakie.vlines!(ax2,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth)

            mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,tr.H0,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

            if type == :existing
                if n==1
                    ax3 = Axis(plot_mut_hist[1,3],title = L"t>H_{0}")
                else
                    ax3 = Axis(plot_mut_hist[1,3])
                end
            else
                ax3 = Axis(plot_mut_hist[2,3])
            end

            min_t = minimum(mut_size)
            max_t = maximum(mut_size)

            norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];

            if type == :existing
                CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[1])
            else
                CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[2])
            end

            CairoMakie.lines!(ax3,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
            CairoMakie.vlines!(ax3,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth)

            # linkxaxes!([ax1,ax2,ax3]...)
            # linkyaxes!([ax1,ax2,ax3]...)

            hidedecorations!(ax1)
            hidedecorations!(ax2)
            hidedecorations!(ax3)
        end

        # colgap!(plot_geno,10)
        # rowgap!(plot_geno,10)
    
        # colgap!(plot_wait, 10)
        # rowgap!(plot_wait, 10)

        # colgap!(plot_mut_hist, 10)
        # rowgap!(plot_mut_hist, 2)

        # colgap!(plot_epi_types, 10)
        # rowgap!(plot_epi_types, 10)

        colgap!(plot_geno,Relative(0.01))
        rowgap!(plot_geno,Relative(0.05))
    
        colgap!(plot_wait, Relative(0.01))
        rowgap!(plot_wait, Relative(0.05))

        colgap!(plot_mut_hist, Relative(0.01))
        rowgap!(plot_mut_hist, Relative(0.05))

        colgap!(plot_epi_types, Relative(0.01))
        rowgap!(plot_epi_types, Relative(0.05))
    end

    # labels = [L"\text{RSE}",L"\text{Sign epistasis}",L"\text{No epistasis}",L"\text{Single mutation}"]

    # Legend(fig[top_n+1, 9:13], [PolyElement(color=c) for c in colors], labels, framevisible=false,nbanks = 2,orientation = :horizontal)

    # colors = palette(:viridis, 3)

    # Legend(fig[top_n+1,2:8],
    #     vcat([[wt_s_list[1], wt_l_list[1]]],[PolyElement(color=c) for c in colors]),
    #     vcat([L"\mathbb{E}[\text{time to accept}]"],[L"\text{weight edits : existing}",L"\text{weight edits : new}",L"\text{weight edits : remove}"]),framevisible=false,nbanks = 2,orientation = :horizontal)

    # labels_wait =  vcat([L"\mathbb{E}[\text{time to accept}]"],[L"\text{weight edits : existing}",L"\text{weight edits : new}",L"\text{weight edits : remove}"])

    # labels_wait =  [L"\mathbb{E}[\text{time to accept}]"]

    if wait_time_summary == :mea
        labels_wait =  [L"\mathbb{E}[\text{total generations}]"]
    else
        labels_wait =  [L"\text{total generations - [25%,50%,75%] quantiles}"]
    end

    labels_mut =  [L"\text{weight edits : existing}",L"\text{weight edits : new}"]

    labels_epi  = [L"\text{RSE}",L"\text{Sign epistasis}",L"\text{No epistasis}",L"\text{Single mutation}"]

    labels = vcat(labels_wait,labels_epi)

    symbol_wait = [[wt_s_list[1], wt_l_list[1]]]

    symbol_mut = [PolyElement(color=c) for c in palette(:viridis, 3)[1:2]]

    symbol_epi = [PolyElement(color=c) for c in evo_config.pie_colors]

    symbol_all = vcat(symbol_wait,symbol_epi)

    # Legend(fig[top_n+1, :], symbol_all, labels, framevisible=false,nbanks = 1,orientation = :horizontal,patchsize = (10, 10), rowgap = 10,colgap = 10)

    # Legend(fig[top_n, 4:13, Bottom()], symbol_all, labels, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = 10,colgap = 10)

    legend_row_gap = 2

    Legend(fig[top_n, 2:4, Bottom()], symbol_wait, labels_wait, framevisible=false,nbanks =1,orientation = :horizontal,patchsize = (10, 10),rowgap = legend_row_gap,colgap = 2,padding=(10.0f0, 10.0f0, 0f0, evo_config.fontsize+1.5*legend_row_gap))

    Legend(fig[top_n, 5:9, Bottom()], symbol_mut, labels_mut, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = legend_row_gap,colgap = 2)
    Legend(fig[top_n,  10:13, Bottom()], symbol_epi, labels_epi, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = legend_row_gap,colgap = 2)

    linkyaxes!(ax_wait_list...)
    linkyaxes!(ax_wait_2_list...)

    rowgap!(fig.layout, Relative(0.01))
    colgap!(fig.layout, Relative(0.01))

end

function create_weight_edit_summary!(fig,n,trajectories,mutation_op::MutationOperator,sorted_uep, vertex_top_map, draw_config, node_colors,fontsize,color_scheme)

    weight_names_latex = reshape([L"W_{aa}",L"W_{ab}",L"W_{ac}",L"W_{ba}",L"W_{bb}",L"W_{bc}",L"W_{ca}",L"W_{cb}",L"W_{cc}",L"W_{ma}",L"W_{mb}",L"W_{mc}"],(3,4));

    grid_values = Tuple.(findall(ones(3,4) .> 0))

    colors = reverse(palette(:tab10)[1:4])

    ax_wait_list = []

    plot_geno = fig[grid_values[11]...] = GridLayout()

    ax_geno = Axis(plot_geno[1,1],backgroundcolor = (color_scheme[n],1.),title =L"\text{Minimal Stripe Topology}",aspect = DataAspect())

    draw_grn!(ax_geno,vertex_top_map[sorted_uep[n]],draw_config,node_colors,fontsize,weight_names_latex,true,false)

    norm_type = :pdf

    ###############################

    plot_mut_hist = fig[grid_values[12]...] = GridLayout()

    bins = 50

    min_t_list = [[],[]]
    max_t_list = [[],[]]

    weight_grid_layouts = []

    for (nt,type) in enumerate([:new,:existing])

        if type == :existing
            mut_noise_dist = mutation_op.noise_distribution;
        else
            mut_noise_dist = Uniform(-mutation_op.max_w,mutation_op.max_w);
        end

        mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

        if type == :existing
            if n==1
                ax1 = Axis(plot_mut_hist[1,1],title = L"t<H_{0}")
            else
                ax1 = Axis(plot_mut_hist[1,1])
            end
        else
            ax1 = Axis(plot_mut_hist[2,1])
        end

        if type == :existing
            CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[1])
        else
            CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[2])
        end

        min_t = minimum(mut_size)
        max_t = maximum(mut_size)

        push!(min_t_list[nt],min_t)
        push!(max_t_list[nt],max_t)

        norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];

        CairoMakie.lines!(ax1,LinRange(min_t,max_t,100),norm_pdf,color = :red)
        CairoMakie.vlines!(ax1,0,color = :red,linestyle = "--")

        mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

        if type == :existing
            if n==1
                ax2 = Axis(plot_mut_hist[1,2],title = L"t=H_{0}")
            else
                ax2 = Axis(plot_mut_hist[1,2])
            end
        else
            ax2 = Axis(plot_mut_hist[2,2])
        end

        min_t = minimum(mut_size)
        max_t = maximum(mut_size)

        push!(min_t_list[nt],min_t)
        push!(max_t_list[nt],max_t)

        norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];

        if type == :existing
            CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[1])
        else
            CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[2])
        end

        CairoMakie.lines!(ax2,LinRange(min_t,max_t,100),norm_pdf,color = :red)
        CairoMakie.vlines!(ax2,0,color = :red,linestyle = "--")

        mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,tr.H0+1,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

        if type == :existing
            if n==1
                ax3 = Axis(plot_mut_hist[1,3],title = L"t>H_{0}")
            else
                ax3 = Axis(plot_mut_hist[1,3])
            end
        else
            ax3 = Axis(plot_mut_hist[2,3])
        end

        min_t = minimum(mut_size)
        max_t = maximum(mut_size)

        push!(min_t_list[nt],min_t)
        push!(max_t_list[nt],max_t)

        norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];

        if type == :existing
            CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[1])
        else
            CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[2])
        end

        CairoMakie.lines!(ax3,LinRange(min_t,max_t,100),norm_pdf,color = :red)
        CairoMakie.vlines!(ax3,0,color = :red,linestyle = "--")

        # linkxaxes!([ax1,ax2,ax3]...)
        # linkyaxes!([ax1,ax2,ax3]...)

        hidedecorations!(ax1)
        hidedecorations!(ax2)
        hidedecorations!(ax3)
    end

    colgap!(plot_mut_hist, 10)
    rowgap!(plot_mut_hist, 10)

    push!(weight_grid_layouts,plot_mut_hist)

    ###############################\

    bins = 15

    for weight_id in 1:10

        plot_mut_hist = fig[grid_values[weight_id]...] = GridLayout()

        for (nt,type) in enumerate([:new,:existing])

            if type == :existing
                mut_noise_dist = Normal(0.0,noise_cv);
            else
                mut_noise_dist = Uniform(-max_w,max_w);
            end

            mut_size = reduce(vcat,map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

            if type == :existing
                ax1 = Axis(plot_mut_hist[1,1],title = L"t<H_{0}")
            else
                ax1 = Axis(plot_mut_hist[2,1])
            end

            min_t = min_t_list[nt][1]
            max_t = max_t_list[nt][1]

            if length(mut_size) > 1
                if type == :existing
                    CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[1])
                else
                    CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[2])
                end
            end

            norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];
            CairoMakie.lines!(ax1,LinRange(min_t,max_t,100),norm_pdf,color = :red)
            CairoMakie.vlines!(ax1,0,color = :red,linestyle = "--")

            mut_size = reduce(vcat,map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

            if type == :existing
                ax2 = Axis(plot_mut_hist[1,2],title = L"t=H_{0}")
            else
                ax2 = Axis(plot_mut_hist[2,2])
            end

            min_t = min_t_list[nt][2]
            max_t = max_t_list[nt][2]

            if length(mut_size) > 1
                if type == :existing
                    CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[1])
                else
                    CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[2])
                end
            end

            norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];
            CairoMakie.lines!(ax2,LinRange(min_t,max_t,100),norm_pdf,color = :red)
            CairoMakie.vlines!(ax2,0,color = :red,linestyle = "--")

            mut_size = reduce(vcat,map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,tr.H0+1,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

            if type == :existing
                ax3 = Axis(plot_mut_hist[1,3],title = L"t>H_{0}")
            else
                ax3 = Axis(plot_mut_hist[2,3])
            end

            min_t = min_t_list[nt][3]
            max_t = max_t_list[nt][3]

            if length(mut_size) > 1
                if type == :existing
                    CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[1])
                else
                    CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[2])
                end
            end

            norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];
            CairoMakie.lines!(ax3,LinRange(min_t,max_t,100),norm_pdf,color = :red)
            CairoMakie.vlines!(ax3,0,color = :red,linestyle = "--")

            # linkxaxes!([ax1,ax2,ax3]...)
            # linkyaxes!([ax1,ax2,ax3]...)

            hidedecorations!(ax1)
            hidedecorations!(ax2)
            hidedecorations!(ax3)
        end

        colgap!(plot_mut_hist, 10)
        rowgap!(plot_mut_hist, 10)

        push!(weight_grid_layouts,plot_mut_hist)

    end

    for (label, layout) in zip(weight_names_latex, weight_grid_layouts[2:end])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = 26,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right)
    end

    Label(weight_grid_layouts[1][1, 1, TopLeft()], L"\text{All}",
    fontsize = 26,
    font = :bold,
    padding = (0, 5, 5, 0),
    halign = :right)

    colors = palette(:viridis, 3)

end

function create_weight_edit_summary!(fig,n,trajectories,mutation_op::MutationOperatorDual,sorted_uep, vertex_top_map, draw_config, node_colors,fontsize,color_scheme)

    weight_names_latex = reshape([L"W_{aa}",L"W_{ab}",L"W_{ac}",L"W_{ba}",L"W_{bb}",L"W_{bc}",L"W_{ca}",L"W_{cb}",L"W_{cc}",L"W_{ma}",L"W_{mb}",L"W_{mc}"],(3,4));

    grid_values = Tuple.(findall(ones(3,4) .> 0))

    colors = reverse(palette(:tab10)[1:4])

    ax_wait_list = []

    plot_geno = fig[grid_values[11]...] = GridLayout()

    ax_geno = Axis(plot_geno[1,1],backgroundcolor = (color_scheme[n],1.),title =L"\text{Minimal Stripe Topology}",aspect = DataAspect())

    draw_grn!(ax_geno,vertex_top_map[sorted_uep[n]],draw_config,node_colors,fontsize,weight_names_latex,true,false)

    norm_type = :pdf

    ###############################

    plot_mut_hist = fig[grid_values[12]...] = GridLayout()

    bins = 50

    min_t_list = [[],[]]
    max_t_list = [[],[]]

    weight_grid_layouts = []

    for (nt,type) in enumerate([:new,:existing])

        if type == :existing
            mut_noise_dist = mutation_op.mult_noise_distribution;
        else
            mut_noise_dist = mutation_op.additive_noise_distribution;
        end

        mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

        if type == :existing
            if n==1
                ax1 = Axis(plot_mut_hist[1,1],title = L"t<H_{0}")
            else
                ax1 = Axis(plot_mut_hist[1,1])
            end
        else
            ax1 = Axis(plot_mut_hist[2,1])
        end

        if type == :existing
            CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[1])
        else
            CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[2])
        end

        min_t = minimum(mut_size)
        max_t = maximum(mut_size)

        push!(min_t_list[nt],min_t)
        push!(max_t_list[nt],max_t)

        norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];

        CairoMakie.lines!(ax1,LinRange(min_t,max_t,100),norm_pdf,color = :red)
        CairoMakie.vlines!(ax1,0,color = :red,linestyle = "--")

        mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

        if type == :existing
            if n==1
                ax2 = Axis(plot_mut_hist[1,2],title = L"t=H_{0}")
            else
                ax2 = Axis(plot_mut_hist[1,2])
            end
        else
            ax2 = Axis(plot_mut_hist[2,2])
        end

        min_t = minimum(mut_size)
        max_t = maximum(mut_size)

        push!(min_t_list[nt],min_t)
        push!(max_t_list[nt],max_t)

        norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];

        if type == :existing
            CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[1])
        else
            CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[2])
        end

        CairoMakie.lines!(ax2,LinRange(min_t,max_t,100),norm_pdf,color = :red)
        CairoMakie.vlines!(ax2,0,color = :red,linestyle = "--")

        mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,tr.H0+1,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

        if type == :existing
            if n==1
                ax3 = Axis(plot_mut_hist[1,3],title = L"t>H_{0}")
            else
                ax3 = Axis(plot_mut_hist[1,3])
            end
        else
            ax3 = Axis(plot_mut_hist[2,3])
        end

        min_t = minimum(mut_size)
        max_t = maximum(mut_size)

        push!(min_t_list[nt],min_t)
        push!(max_t_list[nt],max_t)

        norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];

        if type == :existing
            CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[1])
        else
            CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[2])
        end

        CairoMakie.lines!(ax3,LinRange(min_t,max_t,100),norm_pdf,color = :red)
        CairoMakie.vlines!(ax3,0,color = :red,linestyle = "--")

        # linkxaxes!([ax1,ax2,ax3]...)
        # linkyaxes!([ax1,ax2,ax3]...)

        hidedecorations!(ax1)
        hidedecorations!(ax2)
        hidedecorations!(ax3)
    end

    colgap!(plot_mut_hist, 10)
    rowgap!(plot_mut_hist, 10)

    push!(weight_grid_layouts,plot_mut_hist)

    ###############################\

    bins = 15

    for weight_id in 1:10

        plot_mut_hist = fig[grid_values[weight_id]...] = GridLayout()

        for (nt,type) in enumerate([:new,:existing])

            if type == :existing
                mut_noise_dist = mutation_op.mult_noise_distribution;
            else
                mut_noise_dist = mutation_op.additive_noise_distribution;
            end    

            mut_size = reduce(vcat,map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

            if type == :existing
                ax1 = Axis(plot_mut_hist[1,1],title = L"t<H_{0}")
            else
                ax1 = Axis(plot_mut_hist[2,1])
            end

            min_t = min_t_list[nt][1]
            max_t = max_t_list[nt][1]

            if length(mut_size) > 1
                if type == :existing
                    CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[1])
                else
                    CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[2])
                end
            end

            norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];
            CairoMakie.lines!(ax1,LinRange(min_t,max_t,100),norm_pdf,color = :red)
            CairoMakie.vlines!(ax1,0,color = :red,linestyle = "--")

            mut_size = reduce(vcat,map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

            if type == :existing
                ax2 = Axis(plot_mut_hist[1,2],title = L"t=H_{0}")
            else
                ax2 = Axis(plot_mut_hist[2,2])
            end

            min_t = min_t_list[nt][2]
            max_t = max_t_list[nt][2]

            if length(mut_size) > 1
                if type == :existing
                    CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[1])
                else
                    CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[2])
                end
            end

            norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];
            CairoMakie.lines!(ax2,LinRange(min_t,max_t,100),norm_pdf,color = :red)
            CairoMakie.vlines!(ax2,0,color = :red,linestyle = "--")

            mut_size = reduce(vcat,map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,tr.H0+1,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

            if type == :existing
                ax3 = Axis(plot_mut_hist[1,3],title = L"t>H_{0}")
            else
                ax3 = Axis(plot_mut_hist[2,3])
            end

            min_t = min_t_list[nt][3]
            max_t = max_t_list[nt][3]

            if length(mut_size) > 1
                if type == :existing
                    CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[1])
                else
                    CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[2])
                end
            end

            norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];
            CairoMakie.lines!(ax3,LinRange(min_t,max_t,100),norm_pdf,color = :red)
            CairoMakie.vlines!(ax3,0,color = :red,linestyle = "--")

            # linkxaxes!([ax1,ax2,ax3]...)
            # linkyaxes!([ax1,ax2,ax3]...)

            hidedecorations!(ax1)
            hidedecorations!(ax2)
            hidedecorations!(ax3)
        end

        colgap!(plot_mut_hist, 10)
        rowgap!(plot_mut_hist, 10)

        push!(weight_grid_layouts,plot_mut_hist)

    end

    for (label, layout) in zip(weight_names_latex, weight_grid_layouts[2:end])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = 26,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right)
    end

    Label(weight_grid_layouts[1][1, 1, TopLeft()], L"\text{All}",
    fontsize = 26,
    font = :bold,
    padding = (0, 5, 5, 0),
    halign = :right)

    colors = palette(:viridis, 3)

end

binomialp(k,p,n) = (factorial(n)/(factorial(k)*factorial(n-k)))*(p^k)*((1-p)^(n-k))

function prob_k_mutations(k,p,n)

    binomialp(k,p,n) / (1 - binomialp(0,p,n))

end

function create_mutation_number_summary!(fig,trajectories,top_n,mut_prob,sorted_uep, vertex_top_map, draw_config, node_colors,fontsize,color_scheme)

    null_mut_n = [prob_k_mutations(k,mut_prob,10) for k in 1:10]

    ax_list = []

    for n in 1:top_n

        plot_geno = fig[n, 1] = GridLayout()
        plot_mut_count = fig[n, 2:6] = GridLayout()

        if n==1
            ax_geno = Axis(plot_geno[1,1],backgroundcolor = (color_scheme[n],1.),title =L"\text{Minimal Stripe Topology}",aspect = DataAspect())
        else
            ax_geno = Axis(plot_geno[1,1],backgroundcolor = (color_scheme[n],1.),aspect = DataAspect())
        end

        # top = Int.(reshape(vertex_top_map[sorted_uep[n]],(3,4)))

        # draw_grn_layout!(ax_geno,top,e_width,vertex_size,arrow_size,arrow_shift,sw,fixed_layout,selfedge_size,node_colors,false)

        draw_grn!(ax_geno,vertex_top_map[sorted_uep[n]],draw_config,node_colors,fontsize,false,false)

        #######################

        if n == 1
            ax_mc_lH0 = Axis(plot_mut_count[1,1],title = L"t<H_{0}", xlabel = L"\text{Number of weight edits per mutation}", ylabel = L"\text{Frequency}")
        else
            ax_mc_lH0 = Axis(plot_mut_count[1,1],xlabel = L"\text{Number of weight edits per mutation}", ylabel = L"\text{Frequency}")
        end

        mut_count_prop = sum(reduce(hcat,map(tr->calculate_mut_type_count(get_mut_n(tr,1,tr.H0-2), [i for i in 1:10]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n]) & (tr.H0-2 > 0),trajectories))),dims = 2)[:,1]

        CairoMakie.barplot!(ax_mc_lH0,[i for i in 1:10],mut_count_prop ./ sum(mut_count_prop))
        CairoMakie.lines!(ax_mc_lH0,[i for i in 1:10],null_mut_n,color = :red)
        CairoMakie.scatter!(ax_mc_lH0,[i for i in 1:10],null_mut_n,color = :red, marker = 'x')

        if n == 1
            ax_mc_H0 = Axis(plot_mut_count[1,2],title = L"t=H_{0}",xlabel = L"\text{Number of weight edits per mutation}", ylabel = L"\text{Frequency}")
        else
            ax_mc_H0 = Axis(plot_mut_count[1,2],xlabel = L"\text{Number of weight edits per mutation}", ylabel = L"\text{Frequency}")
        end

        mut_count_prop = sum(reduce(hcat,map(tr->calculate_mut_type_count(get_mut_n(tr,tr.H0-1,tr.H0-1), [i for i in 1:10]),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n] ,trajectories))),dims = 2)[:,1]

        CairoMakie.barplot!(ax_mc_H0,[i for i in 1:10],mut_count_prop ./ sum(mut_count_prop))
        CairoMakie.lines!(ax_mc_H0,[i for i in 1:10],null_mut_n,color = :red)
        CairoMakie.scatter!(ax_mc_H0,[i for i in 1:10],null_mut_n,color = :red, marker = 'x')

        if n == 1
            ax_mc_uH0 = Axis(plot_mut_count[1,3],title = L"t>H_{0}",xlabel = L"\text{Number of weight edits per mutation}", ylabel = L"\text{Frequency}")
        else
            ax_mc_uH0 = Axis(plot_mut_count[1,3],xlabel = L"\text{Number of weight edits per mutation}", ylabel = L"\text{Frequency}")
        end

        mut_count_prop = sum(reduce(hcat,map(tr->calculate_mut_type_count(get_mut_n(tr,tr.H0+1,length(tr.geno_traj)-1),[i for i in 1:10]),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))),dims = 2)[:,1]

        CairoMakie.barplot!(ax_mc_uH0,[i for i in 1:10],mut_count_prop ./ sum(mut_count_prop))
        CairoMakie.lines!(ax_mc_uH0,[i for i in 1:10],null_mut_n,color = :red)
        CairoMakie.scatter!(ax_mc_uH0,[i for i in 1:10],null_mut_n,color = :red, marker = 'x')

        ax_mc_lH0.xticks = (1:10,string.(1:10))
        ax_mc_H0.xticks = (1:10,string.(1:10))
        ax_mc_uH0.xticks = (1:10,string.(1:10))

        colgap!(plot_mut_count, 10)
        rowgap!(plot_mut_count, 10)

        push!(ax_list,ax_mc_lH0)
        push!(ax_list,ax_mc_H0)
        push!(ax_list,ax_mc_uH0)
    end

    linkyaxes!(ax_list...)
end

function plot_mst_explanation!(fig,trajectories,example_id,draw_config,draw_config_generator,node_colors,fontsize,annotate,twin_axis)

    ax1 = Axis(fig[1,2])
    ax2 = Axis(fig[1,2])
    ax_geno = Axis(fig[1,1],backgroundcolor = RGBf(0.98, 0.98, 0.98))

    development = DefaultGRNSolver()

    mst = trajectories[example_id].minimal_stripe_subgraphs[end]

    orig_network = trajectories[example_id].geno_traj[end]
    mst_network = abs.(mst) .* orig_network

    orig_pheno = Individual(reshape(orig_network,(3,4)),grn_parameters,development).phenotype.u[end]
    mst_pheno = Individual(reshape(mst_network,(3,4)),grn_parameters,development).phenotype.u[end]

    if twin_axis
        CairoMakie.lines!(ax1,orig_pheno[3,:], color = :red, linewidth = 4.)
        CairoMakie.lines!(ax2,mst_pheno[3,:], color = :black, linewidth = 4.)
    else
        CairoMakie.lines!(ax1,orig_pheno[3,:], color = :red, linewidth = 4.)
        CairoMakie.lines!(ax1,mst_pheno[3,:], color = :black, linewidth = 4.)
    end

    draw_grn_mutant!(ax_geno,mst_network,orig_network,draw_config,draw_config_generator,node_colors,fontsize,annotate)

    # draw_grn_layout_mutant!(ax_geno,reshape(orig_network,(3,4)),reshape(mst_network,(3,4)),1.2*e_width,1.5*vertex_size,1.5*arrow_size,arrow_shift,sw,fixed_layout,selfedge_size,node_colors)

    hidedecorations!(ax)

    hidedecorations!(ax1)
    hidedecorations!(ax2)

end

mutable struct pred_config
    fontsize
    perf_linewidth
    entropy_markersize
    entropy_linewidth

    entropy_hist_scale 
end

function create_prediction_summary!(fig,trajectories_p,pred_type,max_ce,pred_config)

    cp = palette(:viridis, 3)

    performance_plot = fig[1, 1:2] = GridLayout()
    entropy_plot = fig[2, 1:2] = GridLayout()

    ax_train = Axis(performance_plot[1,1], xlabel = L"v \text{, cumulative weight edits}", title = L"\text{Train (in-sample)}", ylabel = L"\text{% of sample}")
    ax_test = Axis(performance_plot[1,2],xlabel = L"v \text{, cumulative weight edits}", title = L"\text{Train (out-of-sample)}")

    ax_train_accuracy = Axis(performance_plot[1,1])
    ax_test_accuracy = Axis(performance_plot[1,2],yticklabelcolor = :green, yaxisposition = :right,ylabelcolor = :green)

    hideydecorations!(ax_test)
    hideydecorations!(ax_train_accuracy)

    ax_entropy = Axis(entropy_plot[1,1:2],xlabel = L"v \text{, cumulative weight edits}", ylabel = L"\text{Prediction entropies}")

    all_bar_counts = []
    all_bar_stack = []
    all_bar_x = []

    total_tr = 0
    total_te = 0

    if isnothing(max_ce)
        max_ce = maximum(map(tr->tr.weight_edits[tr.H0],trajectories_p))
    end

    pop = filter(tr->tr.train_test_indicator == :train,trajectories_p)

    mean_accuracies = []
    mean_null_accuracies = []

    for r in 1:max_ce
        
        pop_achieved = filter(tr->v_restricted_label_inclusion(tr,x->weight_edit_restriction_measure(x,r),:H0),pop)

        pop_not_achieved = filter(tr->!v_restricted_label_inclusion(tr,x->weight_edit_restriction_measure(x,r),:H0),pop)

        accuracies  = map(tr->v_restricted_accuracy(tr,x->weight_edit_restriction_measure(x,r),pred_type),pop_not_achieved)

        null_accuracies  = map(tr->v_restricted_accuracy(tr,x->weight_edit_restriction_measure(x,0),pred_type),pop_not_achieved)

        bar_counts = [count(x->x==0,accuracies),count(x->x==1,accuracies),length(pop_achieved)] ./ length(trajectories_p)

        bar_stack = [1,2,3]

        bar_x = [r,r,r]

        push!(all_bar_counts,bar_counts)
        push!(all_bar_stack,bar_stack)
        push!(all_bar_x,bar_x)

        if length(accuracies) > 1
            push!(mean_accuracies,mean(accuracies))
            push!(mean_null_accuracies,mean(null_accuracies))
        end

    end

    train_ac = CairoMakie.lines!(ax_train_accuracy,Float64.(mean_accuracies),color = :green, linestyle = "--",linewidth = pred_config.perf_linewidth)
    train_ac_null = CairoMakie.lines!(ax_train_accuracy,Float64.(mean_null_accuracies),color = :cyan,linestyle = "--",linewidth = pred_config.perf_linewidth)

    all_bar_counts = reduce(vcat,all_bar_counts)
    all_bar_stack = reduce(vcat,all_bar_stack)
    all_bar_x = reduce(vcat,all_bar_x);

    CairoMakie.barplot!(ax_train,all_bar_x,all_bar_counts,stack = all_bar_stack,color = all_bar_stack)

    ################################

    all_bar_counts = []
    all_bar_stack = []
    all_bar_x = []

    total_tr = 0
    total_te = 0

    pop = filter(tr->tr.train_test_indicator == :test,trajectories_p)

    mean_accuracies = []
    mean_null_accuracies = []

    for r in 1:max_ce
        
        pop_achieved = filter(tr->v_restricted_label_inclusion(tr,x->weight_edit_restriction_measure(x,r),:H0),pop)

        pop_not_achieved = filter(tr->!v_restricted_label_inclusion(tr,x->weight_edit_restriction_measure(x,r),:H0),pop)

        accuracies  = map(tr->v_restricted_accuracy(tr,x->weight_edit_restriction_measure(x,r),pred_type),pop_not_achieved)

        null_accuracies  = map(tr->v_restricted_accuracy(tr,x->weight_edit_restriction_measure(x,0),pred_type),pop_not_achieved)

        # mean_ent = map(tr->v_restricted_entropy(tr,x->weight_edit_restriction_measure(x,r),:gt),pop_not_achieved)

        # if length(mean_ent) > 10
        # hist!(ax_entropy, mean_ent, scale_to= 0.6, offset=r, direction=:x)
        # push!(all_mean_ent,mean(mean_ent))
        # end

        bar_counts = [count(x->x==0,accuracies),count(x->x==1,accuracies),length(pop_achieved)] ./ length(trajectories_p)

        bar_stack = [1,2,3]

        bar_x = [r,r,r]

        push!(all_bar_counts,bar_counts)
        push!(all_bar_stack,bar_stack)
        push!(all_bar_x,bar_x)

        if length(accuracies) > 1
            push!(mean_accuracies,mean(accuracies))
            push!(mean_null_accuracies,mean(null_accuracies))
        end

    end

    test_ac = CairoMakie.lines!(ax_test_accuracy,Float64.(mean_accuracies),color = :blue,linestyle = "--",linewidth = pred_config.perf_linewidth)
    test_ac_null = CairoMakie.lines!(ax_test_accuracy,Float64.(mean_null_accuracies),color = :cyan,linestyle = "--",linewidth = pred_config.perf_linewidth)

    all_bar_counts = reduce(vcat,all_bar_counts)
    all_bar_stack = reduce(vcat,all_bar_stack)
    all_bar_x = reduce(vcat,all_bar_x);

    CairoMakie.barplot!(ax_test,all_bar_x,all_bar_counts,stack = all_bar_stack,color = all_bar_stack)

    ####################

    all_mean_ent = []

    for r in 1:max_ce

        pop_not_achieved = filter(tr->!v_restricted_label_inclusion(tr,x->weight_edit_restriction_measure(x,r),:H0),trajectories_p)
        
        mean_ent = map(tr->v_restricted_entropy(tr,x->weight_edit_restriction_measure(x,r),pred_type),pop_not_achieved)

        hist!(ax_entropy, mean_ent, scale_to=pred_config.entropy_hist_scale, offset=r, direction=:x, bins = 25, color = (:grey,0.5))
        push!(all_mean_ent,mean(mean_ent))

        # 0.6

    end

    ent_line = CairoMakie.lines!(ax_entropy,Float64.(all_mean_ent),color = :red, linewidth = pred_config.entropy_linewidth)
    ent_scatt = CairoMakie.scatter!(ax_entropy,Float64.(all_mean_ent),color = :red, markersize = pred_config.entropy_markersize)

    ##################

    # linkyaxes!(ax_train,ax_test)
    # linkyaxes!(ax_train_accuracy,ax_test_accuracy)

    linkyaxes!([ax_train,ax_test,ax_train_accuracy,ax_test_accuracy,ax_entropy]...)

    valid_xticks = filter(x->x%2==1,1:max_ce)

    ax_train.xticks = (valid_xticks,string.(valid_xticks))
    ax_test.xticks = (valid_xticks,string.(valid_xticks))

    ax_entropy.xticks = (1:max_ce,string.(1:max_ce))

    CairoMakie.hidedecorations!(ax_train,label = false,ticklabels = false,ticks = false,minorticks = false)

    CairoMakie.hidedecorations!(ax_test,label = false,ticklabels = false,ticks = false,minorticks = false)

    hidespines!(ax_train_accuracy)

    hideydecorations!(ax_test_accuracy,label = false,ticklabels = false,ticks = false,minorticks = false)
    hidexdecorations!(ax_test_accuracy)

    hideydecorations!(ax_train_accuracy,label = false,ticklabels = false,ticks = false,minorticks = false)
    hidexdecorations!(ax_train_accuracy)


    Legend(fig[3, 1:2],
        vcat([train_ac,train_ac_null,[ent_line,ent_scatt]],[PolyElement(color=c) for c in cp]),
        vcat([L"v\text{-restricted accuracy}",L"\text{Null accuracy}",L"\text{Average entropy}"],[L"\text{Incorrect}",L"\text{Correct}",L"M_{H_0} \text{ assigned}"]),framevisible=false,orientation = :horizontal,nbanks = 2,
        patchsize = (10, 10),rowgap = 2,colgap = 7,padding=(10.0f0, 10.0f0, 0.0f0, 0.0f0))

    
    colgap!(performance_plot,Relative(0.01))
    rowgap!(performance_plot,Relative(0.05))

    # colgap!(entropy_plot, Relative(0.01))
    # rowgap!(entropy_plot, Relative(0.05))

    rowgap!(fig.layout, Relative(0.05))
    colgap!(fig.layout, Relative(0.01))

end

function create_prediction_summary!(fig,trajectories_p,pred_type,max_ce,predict_label_to_vertex_rev,pred_config)

    cp = palette(:viridis, 3)

    performance_plot = fig[1, 1:2] = GridLayout()
    entropy_plot = fig[2, 1:2] = GridLayout()

    ax_train = Axis(performance_plot[1,1], xlabel = L"v \text{, cumulative weight edits}", title = L"\text{Train (in-sample)}", ylabel = L"\text{% of sample}")
    ax_test = Axis(performance_plot[1,2],xlabel = L"v \text{, cumulative weight edits}", title = L"\text{Train (out-of-sample)}")

    ax_train_accuracy = Axis(performance_plot[1,1])
    ax_test_accuracy = Axis(performance_plot[1,2],yticklabelcolor = :green, yaxisposition = :right,ylabelcolor = :green)

    hideydecorations!(ax_test)
    hideydecorations!(ax_train_accuracy)

    ax_entropy = Axis(entropy_plot[1,1:2],xlabel = L"v \text{, cumulative weight edits}", ylabel = L"\text{Prediction entropies}")

    all_bar_counts = []
    all_bar_stack = []
    all_bar_x = []

    total_tr = 0
    total_te = 0

    if isnothing(max_ce)
        max_ce = maximum(map(tr->tr.weight_edits[tr.H0],trajectories_p))
    end

    pop = filter(tr->tr.train_test_indicator == :train,trajectories_p)

    mean_accuracies = []
    mean_null_accuracies = []
    roc = []

    for r in 1:max_ce
        
        pop_achieved = filter(tr->v_restricted_label_inclusion(tr,x->weight_edit_restriction_measure(x,r),:H0),pop)

        pop_not_achieved = filter(tr->!v_restricted_label_inclusion(tr,x->weight_edit_restriction_measure(x,r),:H0),pop)

        accuracies  = map(tr->v_restricted_accuracy(tr,x->weight_edit_restriction_measure(x,r),pred_type),pop_not_achieved)

        null_accuracies  = map(tr->v_restricted_accuracy(tr,x->weight_edit_restriction_measure(x,0),pred_type),pop_not_achieved)

        pop_na_labels = map(tr->predict_label_to_vertex_rev[tr.inc_metagraph_vertices[tr.H0]],pop_not_achieved)

        pred_prob  = reduce(hcat,map(tr->v_restricted_probabilities(tr,x->weight_edit_restriction_measure(x,r),pred_type),pop_not_achieved)) |> transpose |> collect;

        roc_score = roc_auc_score(pop_na_labels,pred_prob,multi_class = "ovr", average = "weighted")

        bar_counts = [count(x->x==0,accuracies),count(x->x==1,accuracies),length(pop_achieved)] ./ length(trajectories_p)

        bar_stack = [1,2,3]

        bar_x = [r,r,r]

        push!(all_bar_counts,bar_counts)
        push!(all_bar_stack,bar_stack)
        push!(all_bar_x,bar_x)

        if length(accuracies) > 1
            push!(roc,roc_score)
            push!(mean_accuracies,mean(accuracies))
            push!(mean_null_accuracies,mean(null_accuracies))
        end

    end

    train_ac = CairoMakie.lines!(ax_train_accuracy,Float64.(mean_accuracies),color = :green, linestyle = "--",linewidth = pred_config.perf_linewidth)
    train_ac_null = CairoMakie.lines!(ax_train_accuracy,Float64.(mean_null_accuracies),color = :cyan,linestyle = "--",linewidth = pred_config.perf_linewidth)
    train_roc = CairoMakie.lines!(ax_train_accuracy,Float64.(roc),color = :purple,linestyle = "--",linewidth = pred_config.perf_linewidth)

    all_bar_counts = reduce(vcat,all_bar_counts)
    all_bar_stack = reduce(vcat,all_bar_stack)
    all_bar_x = reduce(vcat,all_bar_x);

    CairoMakie.barplot!(ax_train,all_bar_x,all_bar_counts,stack = all_bar_stack,color = all_bar_stack)

    ################################

    all_bar_counts = []
    all_bar_stack = []
    all_bar_x = []

    total_tr = 0
    total_te = 0

    pop = filter(tr->tr.train_test_indicator == :test,trajectories_p)

    mean_accuracies = []
    mean_null_accuracies = []
    roc = []

    for r in 1:max_ce
        
        pop_achieved = filter(tr->v_restricted_label_inclusion(tr,x->weight_edit_restriction_measure(x,r),:H0),pop)

        pop_not_achieved = filter(tr->!v_restricted_label_inclusion(tr,x->weight_edit_restriction_measure(x,r),:H0),pop)

        accuracies  = map(tr->v_restricted_accuracy(tr,x->weight_edit_restriction_measure(x,r),pred_type),pop_not_achieved)

        null_accuracies  = map(tr->v_restricted_accuracy(tr,x->weight_edit_restriction_measure(x,0),pred_type),pop_not_achieved)

        pop_na_labels = map(tr->predict_label_to_vertex_rev[tr.inc_metagraph_vertices[tr.H0]],pop_not_achieved)

        pred_prob  = reduce(hcat,map(tr->v_restricted_probabilities(tr,x->weight_edit_restriction_measure(x,r),pred_type),pop_not_achieved)) |> transpose |> collect;

        roc_score = roc_auc_score(pop_na_labels,pred_prob,multi_class = "ovr", average = "weighted")

        bar_counts = [count(x->x==0,accuracies),count(x->x==1,accuracies),length(pop_achieved)] ./ length(trajectories_p)

        bar_stack = [1,2,3]

        bar_x = [r,r,r]

        push!(all_bar_counts,bar_counts)
        push!(all_bar_stack,bar_stack)
        push!(all_bar_x,bar_x)

        if length(accuracies) > 1
            push!(roc,roc_score)
            push!(mean_accuracies,mean(accuracies))
            push!(mean_null_accuracies,mean(null_accuracies))
        end

    end

    test_ac = CairoMakie.lines!(ax_test_accuracy,Float64.(mean_accuracies),color = :blue,linestyle = "--",linewidth = pred_config.perf_linewidth)
    test_ac_null = CairoMakie.lines!(ax_test_accuracy,Float64.(mean_null_accuracies),color = :cyan,linestyle = "--",linewidth = pred_config.perf_linewidth)
    test_roc = CairoMakie.lines!(ax_test_accuracy,Float64.(roc),color = :purple,linestyle = "--",linewidth = pred_config.perf_linewidth)

    all_bar_counts = reduce(vcat,all_bar_counts)
    all_bar_stack = reduce(vcat,all_bar_stack)
    all_bar_x = reduce(vcat,all_bar_x);

    CairoMakie.barplot!(ax_test,all_bar_x,all_bar_counts,stack = all_bar_stack,color = all_bar_stack)

    ####################

    all_mean_ent = []

    for r in 1:max_ce

        pop_not_achieved = filter(tr->!v_restricted_label_inclusion(tr,x->weight_edit_restriction_measure(x,r),:H0),trajectories_p)
        
        mean_ent = map(tr->v_restricted_entropy(tr,x->weight_edit_restriction_measure(x,r),pred_type),pop_not_achieved)

        hist!(ax_entropy, mean_ent, scale_to=pred_config.entropy_hist_scale, offset=r, direction=:x, bins = 25, color = (:grey,0.5))
        push!(all_mean_ent,mean(mean_ent))

        # 0.6

    end

    ent_line = CairoMakie.lines!(ax_entropy,Float64.(all_mean_ent),color = :red, linewidth = pred_config.entropy_linewidth)
    ent_scatt = CairoMakie.scatter!(ax_entropy,Float64.(all_mean_ent),color = :red, markersize = pred_config.entropy_markersize)

    ##################

    # linkyaxes!(ax_train,ax_test)
    # linkyaxes!(ax_train_accuracy,ax_test_accuracy)

    linkyaxes!([ax_train,ax_test,ax_train_accuracy,ax_test_accuracy,ax_entropy]...)

    valid_xticks = filter(x->x%2==1,1:max_ce)

    ax_train.xticks = (valid_xticks,string.(valid_xticks))
    ax_test.xticks = (valid_xticks,string.(valid_xticks))

    ax_entropy.xticks = (1:max_ce,string.(1:max_ce))

    CairoMakie.hidedecorations!(ax_train,label = false,ticklabels = false,ticks = false,minorticks = false)

    CairoMakie.hidedecorations!(ax_test,label = false,ticklabels = false,ticks = false,minorticks = false)

    hidespines!(ax_train_accuracy)

    hideydecorations!(ax_test_accuracy,label = false,ticklabels = false,ticks = false,minorticks = false)
    hidexdecorations!(ax_test_accuracy)

    hideydecorations!(ax_train_accuracy,label = false,ticklabels = false,ticks = false,minorticks = false)
    hidexdecorations!(ax_train_accuracy)


    Legend(fig[3, 1:2],
        vcat([train_ac,train_ac_null,[ent_line,ent_scatt]],[PolyElement(color=c) for c in cp]),
        vcat([L"v\text{-restricted accuracy}",L"\text{Null accuracy}",L"\text{Average entropy}"],[L"\text{Incorrect}",L"\text{Correct}",L"M_{H_0} \text{ assigned}"]),framevisible=false,orientation = :horizontal,nbanks = 2,
        patchsize = (10, 10),rowgap = 2,colgap = 7,padding=(10.0f0, 10.0f0, 0.0f0, 0.0f0))

    
    colgap!(performance_plot,Relative(0.01))
    rowgap!(performance_plot,Relative(0.05))

    # colgap!(entropy_plot, Relative(0.01))
    # rowgap!(entropy_plot, Relative(0.05))

    rowgap!(fig.layout, Relative(0.05))
    colgap!(fig.layout, Relative(0.01))

end


function create_evo_summary_DM!(fig,trajectories,top_n,mutation_operator,sorted_uep, vertex_top_map,wait_time_summary,evo_config)

    # all_wait_times = reduce(hcat,[cumulative_wait_time(tr) for tr in trajectories]);

    all_wait_times = reduce(hcat,[average_wait_time(tr) for tr in trajectories]);

    ax_wait_list = []

    ax_wait_list = []
    ax_wait_2_list = []

    wt_l_list = []
    wt_s_list = []

    min_t_u = -mutation_operator.max_w
    max_t_u = mutation_operator.max_w

    for n in 1:top_n

        # plot_geno = fig[n, 1] = GridLayout()
        # plot_wait = fig[n, 2] = GridLayout()
        # plot_mut_hist = fig[n, 3:4] = GridLayout()
        # plot_epi_types = fig[n, 5:6] = GridLayout()

        plot_geno = fig[n, 1] = GridLayout()
        plot_wait = fig[n, 2:4] = GridLayout()

        if n==1
            ax_geno = Axis(plot_geno[1,1],backgroundcolor = (evo_config.color_scheme[n],evo_config.color_fade),title =L"M^{(i)}_{N_i}",aspect = DataAspect())
        else
            ax_geno = Axis(plot_geno[1,1],backgroundcolor = (evo_config.color_scheme[n],evo_config.color_fade),aspect = DataAspect())
        end

        # top = Int.(reshape(vertex_top_map[sorted_uep[n]],(3,4)))

        # draw_grn_layout!(ax_geno,top,e_width,vertex_size,arrow_size,arrow_shift,sw,fixed_layout,selfedge_size,node_colors,false)

        draw_grn!(ax_geno,vertex_top_map[sorted_uep[n]],evo_config.draw_config,evo_config.node_colors,evo_config.fontsize,false,false)

        if n == 1
            ax_wait = Axis(plot_wait[1,1], title = L"\mathbb{E}[\text{total weight edits}]",yticklabelsize = 0.8*evo_config.fontsize)
        else
            ax_wait = Axis(plot_wait[1,1],yticklabelsize = 0.8*evo_config.fontsize)
        end

        ax_wait_2 = Axis(plot_wait[1,1], yticklabelcolor = :red, yaxisposition = :right,yscale = log10,yticklabelsize = 0.8*evo_config.fontsize)

        hidespines!(ax_wait_2 )
        hideydecorations!(ax_wait_2,label = false,ticklabels = false,ticks = false,minorticks = false)
        hidexdecorations!(ax_wait_2)

        # #############################

        mut_type_prop_all = []
        mut_type_time_labels = []
        mut_type_labels = []

        # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,1,tr.H0-2), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n]) & (tr.H0-2 > 0),trajectories_p))

        mut_type_prop = map(tr->calculate_mut_type_count(get_mut_type(tr,1,tr.H0-2), [:TF_pm,:TFBS_pm,:TF_dup,:TFBS_dup]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n]) & (tr.H0-2 > 0),trajectories))

        mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

        push!(mut_type_prop_all,mut_type_prop_av)
        push!(mut_type_labels, [1,2,3,4])
        push!(mut_type_time_labels,[1,1,1])

        # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0-1,tr.H0-1), [:existing,:new,:del]),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n] ,trajectories_p))
        mut_type_prop = map(tr->calculate_mut_type_count(get_mut_type(tr,tr.H0-1,tr.H0-1), [:TF_pm,:TFBS_pm,:TF_dup,:TFBS_dup]),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n] ,trajectories))

        mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

        push!(mut_type_prop_all,mut_type_prop_av)
        push!(mut_type_labels, [1,2,3,4])
        push!(mut_type_time_labels,[2,2,2])

        # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0,length(tr.topologies)-1), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n])  & (tr.H0 < length(tr.topologies)),trajectories_p))

        mut_type_prop = map(tr->calculate_mut_type_count(get_mut_type(tr,tr.H0,length(tr.topologies)-1), [:TF_pm,:TFBS_pm,:TF_dup,:TFBS_dup]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n])  & (tr.H0 < length(tr.topologies)),trajectories))

        mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

        push!(mut_type_prop_all,mut_type_prop_av)
        push!(mut_type_labels, [1,2,3,4])
        push!(mut_type_time_labels,[3,3,3])

        mut_type_prop_all = reduce(vcat,mut_type_prop_all)
        mut_type_time_labels = reduce(vcat,mut_type_time_labels)
        mut_type_labels = reduce(vcat,mut_type_labels); 

        CairoMakie.barplot!(ax_wait,mut_type_time_labels,mut_type_prop_all,stack = mut_type_labels,color = mut_type_labels)

        if n == top_n
            ax_wait.xticks = (1:3,[L"t<H_{0}",L"t=H_{0}",L"t>H_{0}" ])
        else
            hidexdecorations!(ax_wait)
        end
        
        #format y ticks to latex numbers

        CairoMakie.hidexdecorations!(ax_wait,label = false,ticklabels = false,ticks = false,minorticks = false)
        CairoMakie.hideydecorations!(ax_wait,label = false,ticklabels = false,ticks = false,minorticks = false,grid = false)

        push!(ax_wait_list,ax_wait)

        ############################

        sample_id = findall(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)

        if wait_time_summary == :mean

            mean_wait = mean(all_wait_times[:,sample_id],dims = 2)[:,1]

            std_error_wait = std(all_wait_times[:,sample_id],dims = 2)[:,1] ./ sqrt(length(sample_id))

            mean_wait_type_labels = [1,2,3]

            wt_l = CairoMakie.lines!(ax_wait_2,mean_wait_type_labels,mean_wait,color = :red,linewidth = evo_config.wait_linewidth)
            wt_s = CairoMakie.scatter!(ax_wait_2,mean_wait_type_labels,mean_wait,color = :red,markersize = evo_config.wait_markersize)

            CairoMakie.errorbars!(ax_wait_2,1:length(mean_wait),mean_wait,5 * std_error_wait,color = :red,whiskerwidth = evo_config.wait_markersize/2)

        else

            median_wait_time = mapslices(row->quantile(row, [0.5]),all_wait_times[:,sample_id],dims =2)[:,1]
            lq_wait_time = mapslices(row->quantile(row, [0.25]),all_wait_times[:,sample_id],dims =2)[:,1]
            uq_wait_time = mapslices(row->quantile(row, [0.75]),all_wait_times[:,sample_id],dims =2)[:,1]

            median_wait_type_labels = [1,2,3]

            wt_l = CairoMakie.lines!(ax_wait_2,median_wait_type_labels,median_wait_time,color = :red,linewidth = evo_config.wait_linewidth)
            wt_s = CairoMakie.scatter!(ax_wait_2,median_wait_type_labels,median_wait_time,color = :red,markersize = evo_config.wait_markersize)

            CairoMakie.rangebars!(ax_wait_2,1:length(median_wait_time),lq_wait_time,uq_wait_time,color = :red,whiskerwidth = evo_config.wait_markersize/2)
        end

        push!(ax_wait_2_list,ax_wait_2)

        push!(wt_l_list,wt_l)
        push!(wt_s_list,wt_s)

        #############################

        colgap!(plot_geno,Relative(0.01))
        rowgap!(plot_geno,Relative(0.05))
    
        colgap!(plot_wait, Relative(0.01))
        rowgap!(plot_wait, Relative(0.05))

        colgap!(plot_mut_hist, Relative(0.01))
        rowgap!(plot_mut_hist, Relative(0.05))

        colgap!(plot_epi_types, Relative(0.01))
        rowgap!(plot_epi_types, Relative(0.05))
    end

    if wait_time_summary == :mea
        labels_wait =  [L"\mathbb{E}[\text{total generations}]"]
    else
        labels_wait =  [L"\text{total generations - [25%,50%,75%] quantiles}"]
    end

    labels_mut =  [L"\text{TF_pm}",L"\text{TFBS_pm}",L"\text{TF_dup}",L"\text{TFBS_dup}"]

    symbol_wait = [[wt_s_list[1], wt_l_list[1]]]

    symbol_mut = [PolyElement(color=c) for c in palette(:viridis, 4)[1:4]]


    # Legend(fig[top_n+1, :], symbol_all, labels, framevisible=false,nbanks = 1,orientation = :horizontal,patchsize = (10, 10), rowgap = 10,colgap = 10)

    # Legend(fig[top_n, 4:13, Bottom()], symbol_all, labels, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = 10,colgap = 10)

    legend_row_gap = 2

    Legend(fig[top_n, 2:4, Bottom()], symbol_wait, labels_wait, framevisible=false,nbanks =1,orientation = :horizontal,patchsize = (10, 10),rowgap = legend_row_gap,colgap = 2,padding=(10.0f0, 10.0f0, 0f0, evo_config.fontsize+1.5*legend_row_gap))

    Legend(fig[top_n, 5:9, Bottom()], symbol_mut, labels_mut, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = legend_row_gap,colgap = 2)

    linkyaxes!(ax_wait_list...)
    linkyaxes!(ax_wait_2_list...)

    rowgap!(fig.layout, Relative(0.01))
    colgap!(fig.layout, Relative(0.01))

end
