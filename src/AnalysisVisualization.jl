# https://juliadatascience.io/makie_layouts

custom_neglog_formatter(values) = map(
    v -> "-10" * Makie.UnicodeFun.to_superscript(round(Int64, v)),
    values
	)

custom_poslog_formatter(values) = map(
    v -> "10" * Makie.UnicodeFun.to_superscript(round(Int64, v)),
    values
    )

custom_log_formatter(values) = map(
    v -> "10" * Makie.UnicodeFun.to_superscript(round(Int64, v)),
    values
    )

function compute_x(x, width, gap, dodge, dodge_gap)
    scale_width(dodge_gap, n_dodge) = (1 - (n_dodge - 1) * dodge_gap) / n_dodge
    function shift_dodge(i, dodge_width, dodge_gap)
        (dodge_width - 1) / 2 + (i - 1) * (dodge_width + dodge_gap)
    end
    width *= 1 - gap
    n_dodge = maximum(dodge)
    dodge_width = scale_width(dodge_gap, n_dodge)
    shifts = shift_dodge.(dodge, dodge_width, dodge_gap)
    return x .+ width .* shifts
end

binomialp(k,p,n) = (factorial(n)/(factorial(k)*factorial(n-k)))*(p^k)*((1-p)^(n-k))

function prob_k_mutations(k,p,n)

    binomialp(k,p,n) / (1 - binomialp(0,p,n))

end

######################

function plot_convergence_rate!(ax,conv_time,n_trials,max_gen)

    cum_conv = [sum(conv_time .< i)/n_trials for i in 1:max_gen];

    CairoMakie.lines!(ax,cum_conv,color = :blue,linewidth = 4.)

    ax.title = "A. Convergence Rate"

    ax.xlabel = "Number of generations"
    ax.ylabel = "% of trajectories converged"

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

    draw_grn_mutant!(ax_geno,mst_network,orig_network,draw_config,draw_config_generator,node_colors,fontsize,:red,annotate,false)

    # draw_grn_layout_mutant!(ax_geno,reshape(orig_network,(3,4)),reshape(mst_network,(3,4)),1.2*e_width,1.5*vertex_size,1.5*arrow_size,arrow_shift,sw,fixed_layout,selfedge_size,node_colors)

    hidedecorations!(ax)

    hidedecorations!(ax1)
    hidedecorations!(ax2)

end

function plot_sf_explanation!(fig,trajectories,example_id,draw_config,draw_config_generator,node_colors,fontsize,annotate,twin_axis)

    ax1 = Axis(fig[1,2])
    ax2 = Axis(fig[1,2])
    ax_geno = Axis(fig[1,1],backgroundcolor = RGBf(0.98, 0.98, 0.98))

    development = DefaultGRNSolver()

    tr = trajectories[example_id]

    no_stripe = tr.geno_traj[tr.H0-1]
    stripe = tr.geno_traj[tr.H0]

    # stripe_network = abs.(stripe) .* no_stripe

    orig_pheno = Individual(reshape(no_stripe,(3,4)),grn_parameters,development).phenotype.u[end]
    stripe_pheno = Individual(reshape(stripe,(3,4)),grn_parameters,development).phenotype.u[end]

    x=LinRange(0,1,Nc)
    if twin_axis
        CairoMakie.lines!(ax1,x,orig_pheno[3,:], color = :black, linewidth = 4.)
        CairoMakie.lines!(ax2,x,stripe_pheno[3,:], color = :red, linewidth = 4., linestyle = :dash)
    else
        CairoMakie.lines!(ax1,x,orig_pheno[3,:], color = :black, linewidth = 4.)
        CairoMakie.lines!(ax1,x,stripe_pheno[3,:], color = :red, linewidth = 4., linestyle = :dash)
    end

    CairoMakie.lines!(ax1,x,y->10*morph(y),color = node_colors[4], linewidth = 4.)

    # draw_grn_mutant!(ax_geno,no_stripe,stripe,draw_config,draw_config_generator,node_colors,fontsize,annotate)

    draw_grn_mutant!(ax_geno,no_stripe,stripe,draw_config,draw_config_generator,node_colors,fontsize,:red,annotate,false)

    # draw_grn_layout_mutant!(ax_geno,reshape(orig_network,(3,4)),reshape(mst_network,(3,4)),1.2*e_width,1.5*vertex_size,1.5*arrow_size,arrow_shift,sw,fixed_layout,selfedge_size,node_colors)

    hidedecorations!(ax_geno)

    hidedecorations!(ax1)
    hidedecorations!(ax2)

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

mutable struct pred_config
    fontsize
    perf_linewidth
    entropy_markersize
    entropy_linewidth

    entropy_hist_scale 
end

function plot_dynamical_summary_portrait!(fig,trajectories,embedding,top_n,minimal_motif_count,sorted_uep,sorted_counts_uep,mst_conf_int,end_parents,vertex_top_map,example_mst,tr_choice,ds_config)

    # trajectories_p_d = filter(tr->tr.inc_metagraph_vertices[end] ∈ sorted_uep[1:top_n],trajectories);

    mo_umap = fig[1:4, 1:4] = GridLayout()

    ex1 = fig[5:7, 1:4] = GridLayout()
    rmh0 = fig[8,1:4] = GridLayout()

    top_n_dict = Dict(v_id=>pos for (pos,v_id) in enumerate(sorted_uep[1:top_n]))

    #### Motif Distribution

    color_sorted_counts_uep = [i <= top_n ? ds_config.color_scheme[i] : :grey for i in 1:length(sorted_counts_uep)]

    # view_sorted_uep_id = sorted_counts_uep .> minimal_motif_count

    view_sorted_uep_id = [i <= top_n for i in 1:length(sorted_counts_uep)]

    # other = mean(sorted_counts_uep[.!view_sorted_uep_id])

    other = sum(sorted_counts_uep[.!view_sorted_uep_id])

    n_norm = sum(sorted_counts_uep)

    sorted_uep_proportions = vcat(sorted_counts_uep[view_sorted_uep_id],[other]) ./ n_norm 

    view_color_sorted_uep = vcat(color_sorted_counts_uep[view_sorted_uep_id],[:grey])

    # conf_int_choices = mst_conf_int[view_sorted_uep_id]

    # push!(conf_int_choices,(minimum(sorted_counts_uep[.!view_sorted_uep_id]) / n_norm,maximum(sorted_counts_uep[.!view_sorted_uep_id]) / n_norm))

    conf_int_choices = copy(mst_conf_int)

    ##############

    # ax1 = Axis(mo_umap[1:2,1:top_n],title = L"\text{Top %$top_n }" * string(top_n) * " MST : " * string(sum(sorted_counts_uep[1:top_n])) * " trajectories", xlabel = L"\text{Dynamics: UMAP 1}", ylabel = L"\text{Dynamics: UMAP 2}")

    count_top_n = round(sum(sorted_uep_proportions[1:top_n])*100, digits = 2)

    ax1 = Axis(mo_umap[3:4,1:top_n], xlabel = L"\text{Dyn: UMAP 1}", ylabel = L"\text{Dyn: UMAP 2}")

    CairoMakie.scatter!(ax1,embedding, color = [haskey(top_n_dict,i) ? (ds_config.color_scheme[top_n_dict[i]],0.5) : (:grey,0.5) for i in end_parents],markersize = ds_config.embed_markersize)

    hidedecorations!(ax1, label = false)

    for i in 1:top_n

        ax_geno = Axis(mo_umap[2,i], backgroundcolor = (ds_config.color_scheme[i],ds_config.color_fade),aspect = DataAspect())

        draw_grn!(ax_geno,vertex_top_map[sorted_uep[i]],ds_config.draw_config,ds_config.node_colors,ds_config.fontsize,false,false)
    end

    #######################

    ax_mo = Axis(mo_umap[1,1:top_n],ylabel  = L"\text{Probabilty}", xlabel = L"M^{(i)}_{N_i}",title = L"\text{Top %$top_n  M^{(i)}_{N_i} : %$count_top_n % of trajectories}")

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

    ax_fitness = Axis(ex1[1:2,1:length(tr_phenotypes)],xlabel = L"\text{Generation}",ylabel = L"\text{Fitness}")

    hideydecorations!(ax_fitness,label = false,ticklabels = false,ticks = false,minorticks = false)

    rline = CairoMakie.lines!(ax_fitness,tr_fd_refine, color = :grey, linewidth = ds_config.fitness_linewidth)
    cline = CairoMakie.lines!(ax_fitness,tr_fd_coarse, linestyle = "--", color = :blue,linewidth = ds_config.fitness_linewidth)

    CairoMakie.scatter!(ax_fitness,tr_top_fitness_rf, color = [progression_cs[i] for i in 1:length(tr_top_fitness_rf)], markersize = ds_config.fitness_markersize, marker = '★')

    h0 = tr_top_fitness_id[minimum(findall(tr_top_stripe_id .== 1))]
    Ni = tr_top_fitness_id[end]

    if h0 != Ni
        v = Int.(floor((h0+Ni)/2))
        ax_fitness.xticks = ([1,h0,v,Ni],[L"1",L"S_0",L"%$v",L"N_i"])
    else
        ax_fitness.xticks = ([1,h0],[L"1",L"S_0 = N_i"])
    end

    # axislegend(ax_fitness, [rline, cline], [L"\mathcal{F}_R", L"\mathcal{F}_C"], position = :rb,
    # orientation = :vertical, labelsize = 10., rowgap = 2., framevisible = false, valign = :bottom)

    Legend(ex1[1:2,1:length(tr_phenotypes)],  [rline, cline], [L"\mathcal{F}_R(\phi)", L"\mathcal{F}_S(\phi)"], framevisible=false,orientation = :vertical,patchsize = (10, 10),rowgap = 2,halign = :right, valign = :bottom)

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

    for n in 1:top_n+1

        if n == top_n+1
            pop = filter(tr->!(tr.inc_metagraph_vertices[end] ∈  sorted_uep[1:top_n]),trajectories)
        else
            pop = filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)
        end

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

    CairoMakie.barplot!(ax_rh0,x,proportions,color = [n==top_n+1 ? :grey : ds_config.color_scheme[n] for n in dodge],dodge = dodge)

    CairoMakie.hidedecorations!(ax_rh0,label = false,ticklabels = false,ticks = false,minorticks = false)

    ax_rh0.xticks = (1:4,[L"M^{(i)}_{S_{0}} = M^{(i)}_{N_i}",L"M^{(i)}_{S_{0}} \subset M^{(i)}_{N_i}",L"M^{(i)}_{N_i} \subset M^{(i)}_{S_{0}}",L"\text{MST change}"])

    CairoMakie.ylims!(ax_rh0,0.,1.)

    for (label, layout) in zip(["A.i", "B.i", "C"], [mo_umap, ex1, rmh0])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = ds_config.caption_fontsize,
            font = :bold,
            padding = (0,ds_config.caption_padding, ds_config.caption_padding, 0),
            halign = :right)
    end

    Label(mo_umap[3, 1, TopLeft()], "A.ii",
    fontsize = ds_config.caption_fontsize,
    font = :bold,
    padding = (0,ds_config.caption_padding, ds_config.caption_padding, 0),
    halign = :right)

    Label(ex1[3, 1, TopLeft()], "B.ii",
    fontsize = ds_config.caption_fontsize,
    font = :bold,
    padding = (0,ds_config.caption_padding, ds_config.caption_padding, 0),
    halign = :right)
    
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

function create_epi_summary_portrait_bar!(fig,trajectories,top_n,mutation_operator::MutationOperatorDual,mut_prob,sorted_uep, vertex_top_map,wait_time_summary,evo_config)

    all_wait_times = reduce(hcat,[average_wait_time(tr) for tr in trajectories]);

    ax_wait_list = []

    ax_wait_list = []
    ax_wait_2_list = []

    wt_l_list = []
    wt_s_list = []

    min_t_u = -mutation_operator.max_w
    max_t_u = mutation_operator.max_w

    # grid_values = Tuple.(findall(ones(Int(floor(sqrt(top_n))),Int(floor(sqrt(top_n)))) .> 0))

    # geno_subplot = fig[1:top_n,1] = GridLayout()
    # wait_subplot = fig[1:top_n,5:7] = GridLayout()
    # epi_subplot = fig[1:top_n,2:4] = GridLayout()

    geno_subplot = fig[1:top_n,1:2] = GridLayout()
    wait_subplot = fig[1:top_n,7:9] = GridLayout()
    epi_subplot = fig[1:top_n,3:6] = GridLayout()

    null_mut_n = [prob_k_mutations(k,mut_prob,10) for k in 1:10];

    ax_epi_list = []

    for n in 1:top_n

        sub_plot = fig[n,:] = GridLayout()

        if n==1
            ax_geno = Axis(geno_subplot[n,1],backgroundcolor = (evo_config.color_scheme[n],evo_config.color_fade),title =L"M^{(i)}_{N_i}",aspect = DataAspect())
        else
            ax_geno = Axis(geno_subplot[n,1],backgroundcolor = (evo_config.color_scheme[n],evo_config.color_fade),aspect = DataAspect())
        end

        draw_grn!(ax_geno,vertex_top_map[sorted_uep[n]],evo_config.draw_config,evo_config.node_colors,evo_config.fontsize,false,false)

        if n ==1
            ax_wait = Axis(wait_subplot[n,:],yticklabelsize = 0.8*evo_config.fontsize,yaxisposition = :right, xticklabelsvisible = false, yticksize= 0.25*evo_config.fontsize, title = L"t<S_{0} \text{  }  t=S_{0}  \text{  } t=S_{0}")
        else
            ax_wait = Axis(wait_subplot[n,:],yticklabelsize = 0.8*evo_config.fontsize,yaxisposition = :right, xticklabelsvisible = false, yticksize= 0.25*evo_config.fontsize)
        end

        ax_wait_2 = Axis(wait_subplot[n,:], yticklabelcolor = :red,yscale = log10,yticklabelsize = 0.8*evo_config.fontsize,yticksize= 0.25*evo_config.fontsize)

        hidespines!(ax_wait_2 )
        hideydecorations!(ax_wait_2,label = false,ticklabels = false,ticks = false,minorticks = false)
        hidexdecorations!(ax_wait_2)

        # #############################

        mut_type_prop_all = []
        mut_type_time_labels = []
        mut_type_labels = []

        # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,1,tr.H0-2), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n]) & (tr.H0-2 > 0),trajectories_p))

        mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,1,tr.H0-2), [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n]) & (tr.H0-2 > 0),trajectories))

        mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

        push!(mut_type_prop_all,mut_type_prop_av)
        push!(mut_type_labels, [1,2,3,4])
        push!(mut_type_time_labels,[1,1,1,1])

        # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0-1,tr.H0-1), [:existing,:new,:del]),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n] ,trajectories_p))
        mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0-1,tr.H0-1), [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n] ,trajectories))

        mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

        push!(mut_type_prop_all,mut_type_prop_av)
        push!(mut_type_labels, [1,2,3,4])
        push!(mut_type_time_labels,[2,2,2,2])

        # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0,length(tr.topologies)-1), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n])  & (tr.H0 < length(tr.topologies)),trajectories_p))

        mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0,length(tr.topologies)-1), [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n])  & (tr.H0 < length(tr.topologies)),trajectories))

        mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

        push!(mut_type_prop_all,mut_type_prop_av)
        push!(mut_type_labels, [1,2,3,4])
        push!(mut_type_time_labels,[3,3,3,3])

        mut_type_prop_all = reduce(vcat,mut_type_prop_all)

        mut_type_time_labels = reduce(vcat,mut_type_time_labels)
        mut_type_labels = reduce(vcat,mut_type_labels); 

        CairoMakie.barplot!(ax_wait,mut_type_time_labels,mut_type_prop_all,stack = mut_type_labels,color = mut_type_labels)

        # ax_wait.xticks = (1:3,[L"t<H_{0}",L"t=H_{0}",L"t>H_{0}" ])

        #format y ticks to latex numbers

        # CairoMakie.hidexdecorations!(ax_wait,label = false,ticklabels = false,ticks = false,minorticks = false)
        CairoMakie.hidexdecorations!(ax_wait)
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

        if n==top_n
            ax_epi = Axis(epi_subplot[n,1],yticklabelsize = 10.,xticklabelsize = 10,xlabel = L"\text{weight changes per mutant}")
        else
            ax_epi = Axis(epi_subplot[n,1],yticklabelsize = 10.,xticklabelsize = 10)
        end

        epi_counts_lS0 = reduce(vcat,map(tr->tr.epistasis[1:tr.H0-2],filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n]) & (tr.H0-2 > 0),trajectories)));
        mut_n_counts_lS0 = reduce(vcat,map(tr->get_mut_n(tr,1,tr.H0-2),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n]) & (tr.H0-2 > 0),trajectories)));

        epi_mutn_counts_lS0 = countmap(zip(epi_counts_lS0,mut_n_counts_lS0) |> collect)
        total_epi_mutn_lS0 = sum(values(epi_mutn_counts_lS0))
        epi_mutn_prop_lS0 = Dict(key=>value/total_epi_mutn_lS0 for (key,value) in epi_mutn_counts_lS0);

        grp_epi_lS0 = reduce(vcat,[[1,2,3,4] for i in 1:10])
        grp_epi_label_lS0 = reduce(vcat,[[:rse,:se,:ne,:sm] for i in 1:10])
        grp_mutn_lS0 = reduce(vcat,[[i,i,i,i] for i in 1:10])

        values_epi_mutn_lS0 = [haskey(epi_mutn_prop_lS0,key) ? epi_mutn_prop_lS0[key] : 0. for key in zip(grp_epi_label_lS0,grp_mutn_lS0)]
        counts_epi_mutn_lS0 = [haskey(epi_mutn_counts_lS0,key) ? epi_mutn_counts_lS0[key] : 0. for key in zip(grp_epi_label_lS0,grp_mutn_lS0)]

        dodge_epi_mutn_lS0 = [1 for _ in values_epi_mutn_lS0]
        color_epi_mutn_lS0 = [evo_config_12.pie_colors[i] for i in grp_epi_lS0];

        #############

        epi_counts_S0 = reduce(vcat,map(tr->tr.epistasis[tr.H0-1],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)));
        mut_n_counts_S0 = reduce(vcat,map(tr->get_mut_n(tr,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)));

        epi_mutn_counts_S0 = countmap(zip(epi_counts_S0,mut_n_counts_S0) |> collect)
        total_epi_mutn_S0 = sum(values(epi_mutn_counts_S0))
        epi_mutn_prop_S0 = Dict(key=>value/total_epi_mutn_S0 for (key,value) in epi_mutn_counts_S0);

        grp_epi_S0 = reduce(vcat,[[1,2,3,4] for i in 1:10])
        grp_epi_label_S0 = reduce(vcat,[[:rse,:se,:ne,:sm] for i in 1:10])
        grp_mutn_S0 = reduce(vcat,[[i,i,i,i] for i in 1:10])

        values_epi_mutn_S0 = [haskey(epi_mutn_prop_S0,key) ? epi_mutn_prop_S0[key] : 0. for key in zip(grp_epi_label_S0,grp_mutn_S0)]
        counts_epi_mutn_S0 = [haskey(epi_mutn_counts_S0,key) ? epi_mutn_counts_S0[key] : 0. for key in zip(grp_epi_label_S0,grp_mutn_S0)]

        dodge_epi_mutn_S0 = [2 for _ in values_epi_mutn_S0]
        color_epi_mutn_S0 = [evo_config_12.pie_colors[i] for i in grp_epi_S0];

        ##############

        epi_counts_hS0 = reduce(vcat,map(tr->tr.epistasis[tr.H0:end],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)));
        mut_n_counts_hS0 = reduce(vcat,map(tr->get_mut_n(tr,tr.H0,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)));

        epi_mutn_counts_hS0 = countmap(zip(epi_counts_hS0,mut_n_counts_hS0) |> collect)
        total_epi_mutn_hS0 = sum(values(epi_mutn_counts_hS0))
        epi_mutn_prop_hS0 = Dict(key=>value/total_epi_mutn_hS0 for (key,value) in epi_mutn_counts_hS0);

        grp_epi_hS0 = reduce(vcat,[[1,2,3,4] for i in 1:10])
        grp_epi_label_hS0 = reduce(vcat,[[:rse,:se,:ne,:sm] for i in 1:10])
        grp_mutn_hS0 = reduce(vcat,[[i,i,i,i] for i in 1:10])

        values_epi_mutn_hS0 = [haskey(epi_mutn_prop_hS0,key) ? epi_mutn_prop_hS0[key] : 0. for key in zip(grp_epi_label_hS0,grp_mutn_hS0)]
        counts_epi_mutn_hS0 = [haskey(epi_mutn_counts_hS0,key) ? epi_mutn_counts_hS0[key] : 0. for key in zip(grp_epi_label_hS0,grp_mutn_hS0)]

        dodge_epi_mutn_hS0 = [3 for _ in values_epi_mutn_hS0]
        color_epi_mutn_hS0 = [evo_config.pie_colors[i] for i in grp_epi_hS0];

        values_epi_mutn_all = reduce(vcat,[values_epi_mutn_lS0,values_epi_mutn_S0,values_epi_mutn_hS0])
        counts_epi_mutn_all = Int.(reduce(vcat,[counts_epi_mutn_lS0,counts_epi_mutn_S0,counts_epi_mutn_hS0]))

        x_epi_mutn_all = reduce(vcat,[grp_mutn_lS0,grp_mutn_S0,grp_mutn_hS0])
        stack_epi_mutn_all = reduce(vcat,[grp_epi_lS0,grp_epi_S0,grp_epi_hS0])
        color_epi_mutn_all = reduce(vcat,[color_epi_mutn_lS0,color_epi_mutn_S0,color_epi_mutn_hS0])
        dodge_epi_mutn_all = reduce(vcat,[dodge_epi_mutn_lS0,dodge_epi_mutn_S0,dodge_epi_mutn_hS0]);

        dodge_gap = 0.03 #Makie default
        gap = 0.2 #Makie default
        width = 1 #Makie default

        xerr = compute_x(x_epi_mutn_all, width, gap, dodge_epi_mutn_all, dodge_gap)[1:4:end];

        values_epi_mutn_totals = [sum(values_epi_mutn_all[4*(n-1)+1:4*n]) for (n,i) in enumerate(1:4:length(values_epi_mutn_all))];
        counts_epi_mutn_totals = [sum(counts_epi_mutn_all[4*(n-1)+1:4*n]) for (n,i) in enumerate(1:4:length(counts_epi_mutn_all))];

        epi_mutn_pvalue = 0.01
        err_values_epi_mutn = reduce(vcat,[confint(MultinomialLRTest(counts_epi_mutn_totals[1:10]),epi_mutn_pvalue),confint(MultinomialLRTest(counts_epi_mutn_totals[11:20]),epi_mutn_pvalue),confint(MultinomialLRTest(counts_epi_mutn_totals[21:30]),epi_mutn_pvalue)]);

        scatter_height_epi_mutn = [maximum([values_epi_mutn_totals[i],values_epi_mutn_totals[10+i],values_epi_mutn_totals[20+i]]) for i in 1:10]

        scatter_height_epi_mutn = reduce(vcat,[scatter_height_epi_mutn for i in 1:3]);
        scatter_epi_mutn_mark = reduce(vcat,[[:dtriangle for i in 1:10],[:star5 for i in 1:10],[:utriangle for i in 1:10]]);
        scatter_epi_mutn_color = reduce(vcat,[[:cyan for i in 1:10],[:orange for i in 1:10],[:purple for i in 1:10]]);

        CairoMakie.barplot!(ax_epi,x_epi_mutn_all,values_epi_mutn_all, color = color_epi_mutn_all, stack = stack_epi_mutn_all, dodge = dodge_epi_mutn_all)

        # rangebars!(ax_epi,xerr,first.(err_values_epi_mutn),last.(err_values_epi_mutn ); whiskerwidth = 8)

        CairoMakie.lines!(ax_epi,null_mut_n, color = :black, linestyle = :dash, linewidth = evo_config.wait_linewidth)
        CairoMakie.scatter!(ax_epi,null_mut_n, color = :black,marker = 'x',markersize = evo_config.wait_markersize)

        CairoMakie.scatter!(ax_epi,xerr, max.(1.4 .* scatter_height_epi_mutn,0.1),color = scatter_epi_mutn_color,marker = scatter_epi_mutn_mark,markersize = evo_config.wait_markersize)

        # max_we = maximum(x_epi_mutn_all[findall(x->x!=0,values_epi_mutn_all)])

        max_we = 5

        CairoMakie.xlims!(ax_epi, 2*minimum(xerr) - 1,max_we + 2*(1-minimum(xerr)))

        # axislegend!()

        ax_epi.xticks = (1:max_we, string.(1:max_we))

        if n != top_n
            hidexdecorations!(ax_epi)
        end

        push!(ax_epi_list,ax_epi)

        ###############################\

    end

    linkyaxes!(ax_epi_list...)

    colgap!(geno_subplot, Relative(0.01))
    rowgap!(geno_subplot, Relative(0.05))

    colgap!(wait_subplot, Relative(0.01))
    rowgap!(wait_subplot, Relative(0.05))

    colgap!(epi_subplot, Relative(0.01))
    rowgap!(epi_subplot, Relative(0.01))

    ylabelwl = Label(epi_subplot[1:top_n,1,Left()], L"\text{Proportion of mutants}", rotation = pi/2, padding = (2.,30.,0.,0.))

    ylabelwl = Label(wait_subplot[1:top_n,1,Left()], L"\text{Average wait time}", rotation = pi/2, padding = (2.,30.,0.,0.), color = :red)
    ylabelwr = Label(wait_subplot[1:top_n,end,Right()], L"\text{Mutant composition}", rotation = pi/2, padding = (30.,0.,0.,0.))

    # title_wt = Label(wait_subplot[1,1:end,TopLeft()], L"\text{Title}", padding = (0.,0.,10.,5.))

    if wait_time_summary == :mean
        labels_wait =  [L"\mathbb{E}[\text{Total generations}]"]
    else
        labels_wait =  [L"\text{Total generations}"]
    end

    labels_mut =  [L"\text{new:+}",L"\text{existing:+}",L"\text{new:} \times",L"\text{existing:} \times"]

    labels_epi  = [L"\text{TD}",L"\text{SD}",L"\text{TI}",L"\text{SIC}"]

    labels_per = [L"t<S_0",L"t=S_0",L"t>S_0"]

    labels = reduce(vcat,[labels_wait,labels_epi,labels_mut])

    symbol_wait = [[wt_s_list[1], wt_l_list[1]]]

    symbol_mut = [PolyElement(color=c) for c in palette(:viridis, 4)[1:4]]

    symbol_epi = [PolyElement(color=c) for c in evo_config.pie_colors]

    symbol_per = [MarkerElement(color = :cyan, marker = :dtriangle, markersize = evo_config.wait_markersize),MarkerElement(color = :orange, marker = :star5, markersize = evo_config.wait_markersize),MarkerElement(color = :purple, marker = :utriangle, markersize = evo_config.wait_markersize)]

    symbol_wait_all = reduce(vcat,[symbol_wait,symbol_mut])
    label_wait_all = reduce(vcat,[labels_wait,labels_mut])

    legend_row_gap = 3

    # Legend(epi_subplot[top_n, :, Bottom()], symbol_epi, labels_epi, framevisible=false,nbanks = 1,orientation = :horizontal,patchsize = (10, 10), colgap = 4, padding=(0.,0.,0f0, evo_config.fontsize+2.5*legend_row_gap))

    Legend(epi_subplot[1, :, Top()], symbol_epi, labels_epi, framevisible=false,nbanks = 1,orientation = :horizontal,patchsize = (10, 10), colgap = 4, padding=(0.,0.,legend_row_gap,0.))

    Legend(wait_subplot[top_n, :, Bottom()], symbol_mut, labels_mut, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), colgap = 4, rowgap = 4, padding=(0.,0.,0f0, evo_config.fontsize+legend_row_gap))

    Legend(geno_subplot[top_n, :, Bottom()], symbol_per, labels_per, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), colgap = 4, rowgap = 4, padding=(0.,0.,0f0, evo_config.fontsize+legend_row_gap))

    linkyaxes!(ax_wait_list...)
    linkyaxes!(ax_wait_2_list...)

    rowgap!(fig.layout, Relative(0.1))
    colgap!(fig.layout, Relative(0.01))

end

function create_epi_single_portrait_bar_v1!(fig,trajectories,mut_prob,all_pheno_mut_chf,wait_time_summary,evo_config)

    all_wait_times = reduce(hcat,[average_wait_time_ext(tr) for tr in trajectories]);

    ax_wait_list = []

    ax_wait_list = []
    ax_wait_2_list = []

    wt_l_list = []
    wt_s_list = []


    # grid_values = Tuple.(findall(ones(Int(floor(sqrt(top_n))),Int(floor(sqrt(top_n)))) .> 0))

    # geno_subplot = fig[1:top_n,1] = GridLayout()
    # wait_subplot = fig[1:top_n,5:7] = GridLayout()
    # epi_subplot = fig[1:top_n,2:4] = GridLayout()

    wait_subplot = fig[1,6:8] = GridLayout()
    epi_subplot = fig[1,1:5] = GridLayout()

    pheno_subplot = fig[2,1:8] = GridLayout()

    null_mut_n = [prob_k_mutations(k,mut_prob,10) for k in 1:10];

    Label

    ax_epi_list = []

    ax_wait = Axis(wait_subplot[1,1],yticklabelsize = 0.8*evo_config.fontsize,yaxisposition = :right, xticklabelsvisible = false, yticksize= 0.25*evo_config.fontsize, title = L"n=1 \text{|} 1<n<S_{0} \text{|}  n=S_{0}  \text{|} n=S_{0}")

    ax_wait_2 = Axis(wait_subplot[1,1], yticklabelcolor = :red,yscale = log10,yticklabelsize = 0.8*evo_config.fontsize,yticksize= 0.25*evo_config.fontsize)

    hidespines!(ax_wait_2 )
    hideydecorations!(ax_wait_2,label = false,ticklabels = false,ticks = false,minorticks = false)
    hidexdecorations!(ax_wait_2)

    # #############################

    mut_type_prop_all = []
    mut_type_time_labels = []
    mut_type_labels = []

    # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,1,tr.H0-2), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n]) & (tr.H0-2 > 0),trajectories_p))

    mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,1,1), [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]),filter(tr->(tr.H0-2 > 0),trajectories))

    mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

    push!(mut_type_prop_all,mut_type_prop_av)
    push!(mut_type_labels, [1,2,3,4])
    push!(mut_type_time_labels,[1,1,1,1])

    # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,1,tr.H0-2), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n]) & (tr.H0-2 > 0),trajectories_p))

    mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,2,tr.H0-2), [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]),filter(tr->(tr.H0-2 > 2),trajectories))

    mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

    push!(mut_type_prop_all,mut_type_prop_av)
    push!(mut_type_labels, [1,2,3,4])
    push!(mut_type_time_labels,[2,2,2,2])

    # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0-1,tr.H0-1), [:existing,:new,:del]),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n] ,trajectories_p))
    mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0-1,tr.H0-1), [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]),trajectories)

    mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

    push!(mut_type_prop_all,mut_type_prop_av)
    push!(mut_type_labels, [1,2,3,4])
    push!(mut_type_time_labels,[3,3,3,3])

    # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0,length(tr.topologies)-1), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n])  & (tr.H0 < length(tr.topologies)),trajectories_p))

    mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0,length(tr.topologies)-1), [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]),filter(tr->(tr.H0 < length(tr.topologies)),trajectories))

    mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

    push!(mut_type_prop_all,mut_type_prop_av)
    push!(mut_type_labels, [1,2,3,4])
    push!(mut_type_time_labels,[4,4,4,4])

    mut_type_prop_all = reduce(vcat,mut_type_prop_all)

    mut_type_time_labels = reduce(vcat,mut_type_time_labels)
    mut_type_labels = reduce(vcat,mut_type_labels); 

    CairoMakie.barplot!(ax_wait,mut_type_time_labels,mut_type_prop_all,stack = mut_type_labels,color = mut_type_labels)

    # ax_wait.xticks = (1:3,[L"t<H_{0}",L"t=H_{0}",L"t>H_{0}" ])

    #format y ticks to latex numbers

    # CairoMakie.hidexdecorations!(ax_wait,label = false,ticklabels = false,ticks = false,minorticks = false)
    CairoMakie.hidexdecorations!(ax_wait)
    CairoMakie.hideydecorations!(ax_wait,label = false,ticklabels = false,ticks = false,minorticks = false,grid = false)

    ############################ noise_distribut

    if wait_time_summary == :mean

        mean_wait = mean(all_wait_times,dims = 2)[:,1]

        std_error_wait = std(all_wait_times,dims = 2)[:,1] ./ sqrt(length(sample_id))

        mean_wait_type_labels = [1,2,3,4]

        wt_l = CairoMakie.lines!(ax_wait_2,mean_wait_type_labels,mean_wait,color = :red,linewidth = evo_config.wait_linewidth)
        wt_s = CairoMakie.scatter!(ax_wait_2,mean_wait_type_labels,mean_wait,color = :red,markersize = evo_config.wait_markersize)

        CairoMakie.errorbars!(ax_wait_2,1:length(mean_wait),mean_wait,5 * std_error_wait,color = :red,whiskerwidth = evo_config.wait_markersize/2)

    else

        median_wait_time = mapslices(row->quantile(row, [0.5]),all_wait_times,dims =2)[:,1]
        lq_wait_time = mapslices(row->quantile(row, [0.25]),all_wait_times,dims =2)[:,1]
        uq_wait_time = mapslices(row->quantile(row, [0.75]),all_wait_times,dims =2)[:,1]

        median_wait_type_labels = [1,2,3,4]

        wt_l = CairoMakie.lines!(ax_wait_2,median_wait_type_labels,median_wait_time,color = :red,linewidth = evo_config.wait_linewidth)
        wt_s = CairoMakie.scatter!(ax_wait_2,median_wait_type_labels,median_wait_time,color = :red,markersize = evo_config.wait_markersize)

        CairoMakie.rangebars!(ax_wait_2,1:length(median_wait_time),lq_wait_time,uq_wait_time,color = :red,whiskerwidth = evo_config.wait_markersize/2)
    end

    push!(wt_l_list,wt_l)
    push!(wt_s_list,wt_s)

    #############################

    ax_epi = Axis(epi_subplot[1,1],yticklabelsize = 10.,xticklabelsize = 10,xlabel = L"\text{weight changes per mutant}",ygridvisible = false,xgridvisible = false)

    epi_counts_1 = reduce(vcat,map(tr->tr.epistasis[1][1],filter(tr->(tr.H0 > 2),trajectories)));
    mut_n_counts_1 = reduce(vcat,map(tr->get_mut_n(tr,1,1),filter(tr->(tr.H0 > 2),trajectories)));

    epi_mutn_counts_1 = countmap(zip(epi_counts_1,mut_n_counts_1) |> collect)
    total_epi_mutn_1 = sum(values(epi_mutn_counts_1))
    epi_mutn_prop_1 = Dict(key=>value/total_epi_mutn_1 for (key,value) in epi_mutn_counts_1);

    grp_epi_1 = reduce(vcat,[[1,2,3,4] for i in 1:10])
    grp_epi_label_1 = reduce(vcat,[[:rse,:se,:ne,:sm] for i in 1:10])
    grp_mutn_1 = reduce(vcat,[[i,i,i,i] for i in 1:10])

    values_epi_mutn_1 = [haskey(epi_mutn_prop_1,key) ? epi_mutn_prop_1[key] : 0. for key in zip(grp_epi_label_1,grp_mutn_1)]
    counts_epi_mutn_1 = [haskey(epi_mutn_counts_1,key) ? epi_mutn_counts_1[key] : 0. for key in zip(grp_epi_label_1,grp_mutn_1)]

    dodge_epi_mutn_1 = [1 for _ in values_epi_mutn_1]
    color_epi_mutn_1 = [evo_config_12.pie_colors[i] for i in grp_epi_1];

    ###################

    epi_counts_lS0 = reduce(vcat,map(tr->first.(tr.epistasis[2:tr.H0-2]),filter(tr->(tr.H0-2 > 0),trajectories)));
    mut_n_counts_lS0 = reduce(vcat,map(tr->get_mut_n(tr,2,tr.H0-2),filter(tr->(tr.H0-2 > 0),trajectories)));

    epi_mutn_counts_lS0 = countmap(zip(epi_counts_lS0,mut_n_counts_lS0) |> collect)
    total_epi_mutn_lS0 = sum(values(epi_mutn_counts_lS0))
    epi_mutn_prop_lS0 = Dict(key=>value/total_epi_mutn_lS0 for (key,value) in epi_mutn_counts_lS0);

    grp_epi_lS0 = reduce(vcat,[[1,2,3,4] for i in 1:10])
    grp_epi_label_lS0 = reduce(vcat,[[:rse,:se,:ne,:sm] for i in 1:10])
    grp_mutn_lS0 = reduce(vcat,[[i,i,i,i] for i in 1:10])

    values_epi_mutn_lS0 = [haskey(epi_mutn_prop_lS0,key) ? epi_mutn_prop_lS0[key] : 0. for key in zip(grp_epi_label_lS0,grp_mutn_lS0)]
    counts_epi_mutn_lS0 = [haskey(epi_mutn_counts_lS0,key) ? epi_mutn_counts_lS0[key] : 0. for key in zip(grp_epi_label_lS0,grp_mutn_lS0)]

    dodge_epi_mutn_lS0 = [2 for _ in values_epi_mutn_lS0]
    color_epi_mutn_lS0 = [evo_config_12.pie_colors[i] for i in grp_epi_lS0];

    #############

    epi_counts_S0 = reduce(vcat,map(tr->tr.epistasis[tr.H0-1][1],trajectories));
    mut_n_counts_S0 = reduce(vcat,map(tr->get_mut_n(tr,tr.H0-1,tr.H0-1),trajectories));

    epi_mutn_counts_S0 = countmap(zip(epi_counts_S0,mut_n_counts_S0) |> collect)
    total_epi_mutn_S0 = sum(values(epi_mutn_counts_S0))
    epi_mutn_prop_S0 = Dict(key=>value/total_epi_mutn_S0 for (key,value) in epi_mutn_counts_S0);

    grp_epi_S0 = reduce(vcat,[[1,2,3,4] for i in 1:10])
    grp_epi_label_S0 = reduce(vcat,[[:rse,:se,:ne,:sm] for i in 1:10])
    grp_mutn_S0 = reduce(vcat,[[i,i,i,i] for i in 1:10])

    values_epi_mutn_S0 = [haskey(epi_mutn_prop_S0,key) ? epi_mutn_prop_S0[key] : 0. for key in zip(grp_epi_label_S0,grp_mutn_S0)]
    counts_epi_mutn_S0 = [haskey(epi_mutn_counts_S0,key) ? epi_mutn_counts_S0[key] : 0. for key in zip(grp_epi_label_S0,grp_mutn_S0)]

    dodge_epi_mutn_S0 = [3 for _ in values_epi_mutn_S0]
    color_epi_mutn_S0 = [evo_config_12.pie_colors[i] for i in grp_epi_S0];

    ##############

    epi_counts_hS0 = reduce(vcat,map(tr->first.(tr.epistasis[tr.H0:end]),trajectories));
    mut_n_counts_hS0 = reduce(vcat,map(tr->get_mut_n(tr,tr.H0,length(tr.geno_traj)-1),trajectories));

    epi_mutn_counts_hS0 = countmap(zip(epi_counts_hS0,mut_n_counts_hS0) |> collect)
    total_epi_mutn_hS0 = sum(values(epi_mutn_counts_hS0))
    epi_mutn_prop_hS0 = Dict(key=>value/total_epi_mutn_hS0 for (key,value) in epi_mutn_counts_hS0);

    grp_epi_hS0 = reduce(vcat,[[1,2,3,4] for i in 1:10])
    grp_epi_label_hS0 = reduce(vcat,[[:rse,:se,:ne,:sm] for i in 1:10])
    grp_mutn_hS0 = reduce(vcat,[[i,i,i,i] for i in 1:10])

    values_epi_mutn_hS0 = [haskey(epi_mutn_prop_hS0,key) ? epi_mutn_prop_hS0[key] : 0. for key in zip(grp_epi_label_hS0,grp_mutn_hS0)]
    counts_epi_mutn_hS0 = [haskey(epi_mutn_counts_hS0,key) ? epi_mutn_counts_hS0[key] : 0. for key in zip(grp_epi_label_hS0,grp_mutn_hS0)]

    dodge_epi_mutn_hS0 = [4 for _ in values_epi_mutn_hS0]
    color_epi_mutn_hS0 = [evo_config.pie_colors[i] for i in grp_epi_hS0];

    values_epi_mutn_all = reduce(vcat,[values_epi_mutn_1,values_epi_mutn_lS0,values_epi_mutn_S0,values_epi_mutn_hS0])
    counts_epi_mutn_all = Int.(reduce(vcat,[counts_epi_mutn_1,counts_epi_mutn_lS0,counts_epi_mutn_S0,counts_epi_mutn_hS0]))

    x_epi_mutn_all = reduce(vcat,[grp_mutn_1,grp_mutn_lS0,grp_mutn_S0,grp_mutn_hS0])
    stack_epi_mutn_all = reduce(vcat,[grp_epi_1,grp_epi_lS0,grp_epi_S0,grp_epi_hS0])
    color_epi_mutn_all = reduce(vcat,[color_epi_mutn_1,color_epi_mutn_lS0,color_epi_mutn_S0,color_epi_mutn_hS0])
    dodge_epi_mutn_all = reduce(vcat,[dodge_epi_mutn_1,dodge_epi_mutn_lS0,dodge_epi_mutn_S0,dodge_epi_mutn_hS0]);

    dodge_gap = 0.03 #Makie default
    gap = 0.2 #Makie default
    width = 1 #Makie default

    xerr = compute_x(x_epi_mutn_all, width, gap, dodge_epi_mutn_all, dodge_gap)[1:4:end];

    values_epi_mutn_totals = [sum(values_epi_mutn_all[4*(n-1)+1:4*n]) for (n,i) in enumerate(1:4:length(values_epi_mutn_all))];
    counts_epi_mutn_totals = [sum(counts_epi_mutn_all[4*(n-1)+1:4*n]) for (n,i) in enumerate(1:4:length(counts_epi_mutn_all))];

    epi_mutn_pvalue = 0.01
    err_values_epi_mutn = reduce(vcat,[confint(MultinomialLRTest(counts_epi_mutn_totals[1:10]),epi_mutn_pvalue),confint(MultinomialLRTest(counts_epi_mutn_totals[11:20]),epi_mutn_pvalue),confint(MultinomialLRTest(counts_epi_mutn_totals[21:30]),epi_mutn_pvalue)]);

    scatter_height_epi_mutn = [maximum([values_epi_mutn_totals[i],values_epi_mutn_totals[10+i],values_epi_mutn_totals[20+i],values_epi_mutn_totals[30+i]]) for i in 1:10]

    scatter_height_epi_mutn = reduce(vcat,[scatter_height_epi_mutn for i in 1:4]);
    scatter_epi_mutn_mark = reduce(vcat,[[:dtriangle for i in 1:10],[:star5 for i in 1:10],[:utriangle for i in 1:10]]);
    scatter_epi_mutn_color = reduce(vcat,[[:cyan for i in 1:10],[:orange for i in 1:10],[:purple for i in 1:10]]);

    scatter_epi_mutn_text = reduce(vcat,[[L"n = 1" for i in 1:10],[L"1 < n < S_0" for i in 1:10],[L"n = S_0" for i in 1:10],[L"n > S_0" for i in 1:10]]);

    CairoMakie.barplot!(ax_epi,x_epi_mutn_all,values_epi_mutn_all, color = color_epi_mutn_all, stack = stack_epi_mutn_all, dodge = dodge_epi_mutn_all,width = width)

    # rangebars!(ax_epi,xerr,first.(err_values_epi_mutn),last.(err_values_epi_mutn ); whiskerwidth = 8)

    CairoMakie.lines!(ax_epi,null_mut_n, color = :black, linestyle = :dash, linewidth = evo_config.wait_linewidth)
    CairoMakie.scatter!(ax_epi,null_mut_n, color = :black,marker = 'x',markersize = evo_config.wait_markersize)

    # CairoMakie.scatter!(ax_epi,xerr, max.(1.4 .* scatter_height_epi_mutn,0.1),color = scatter_epi_mutn_color,marker = scatter_epi_mutn_mark,markersize = evo_config.wait_markersize)

    CairoMakie.text!(ax_epi,xerr, max.(1.05 .* scatter_height_epi_mutn,0.1),text = scatter_epi_mutn_text,rotation = pi/2,fontsize = 10.,align = (:left, :center))

    # max_we = maximum(x_epi_mutn_all[findall(x->x!=0,values_epi_mutn_all)])

    max_we = 4

    CairoMakie.xlims!(ax_epi, 2*minimum(xerr) - 1,max_we + 2*(1-minimum(xerr)))
    CairoMakie.ylims!(ax_epi, 0.,1.)

    # axislegend!()

    ax_epi.xticks = (1:max_we, string.(1:max_we))


    ###############################\

    colgap!(wait_subplot, Relative(0.01))
    rowgap!(wait_subplot, Relative(0.05))

    colgap!(epi_subplot, Relative(0.01))
    rowgap!(epi_subplot, Relative(0.01))

    ylabelwl = Label(epi_subplot[1,1,Left()], L"\text{Proportion of mutants}", rotation = pi/2, padding = (2.,30.,0.,0.))

    ylabelwl = Label(wait_subplot[1,1,Left()], L"\text{Average wait time}", rotation = pi/2, padding = (2.,30.,0.,0.), color = :red)
    ylabelwr = Label(wait_subplot[1,end,Right()], L"\text{Mutant composition}", rotation = pi/2, padding = (30.,0.,0.,0.))

    # title_wt = Label(wait_subplot[1,1:end,TopLeft()], L"\text{Title}", padding = (0.,0.,10.,5.))

    if wait_time_summary == :mean
        labels_wait =  [L"\mathbb{E}[\text{Total generations}]"]
    else
        labels_wait =  [L"\text{Total generations}"]
    end

    labels_mut =  [L"\text{new:+}",L"\text{existing:+}",L"\text{new:} \times",L"\text{existing:} \times"]

    labels_epi  = [L"\text{TD}",L"\text{SD}",L"\text{TI}",L"\text{SIC}"]

    labels_per = [L"t<S_0",L"t=S_0",L"t>S_0"]

    labels = reduce(vcat,[labels_wait,labels_epi,labels_mut])

    symbol_wait = [[wt_s_list[1], wt_l_list[1]]]

    symbol_mut = [PolyElement(color=c) for c in palette(:viridis, 4)[1:4]]

    symbol_epi = [PolyElement(color=c) for c in evo_config.pie_colors]

    symbol_per = [MarkerElement(color = :cyan, marker = :dtriangle, markersize = evo_config.wait_markersize),MarkerElement(color = :orange, marker = :star5, markersize = evo_config.wait_markersize),MarkerElement(color = :purple, marker = :utriangle, markersize = evo_config.wait_markersize)]

    symbol_wait_all = reduce(vcat,[symbol_wait,symbol_mut])
    label_wait_all = reduce(vcat,[labels_wait,labels_mut])

    legend_row_gap = 3

    # Legend(epi_subplot[top_n, :, Bottom()], symbol_epi, labels_epi, framevisible=false,nbanks = 1,orientation = :horizontal,patchsize = (10, 10), colgap = 4, padding=(0.,0.,0f0, evo_config.fontsize+2.5*legend_row_gap))

    Legend(epi_subplot[1, :, Top()], symbol_epi, labels_epi, framevisible=false,nbanks = 1,orientation = :horizontal,patchsize = (10, 10), colgap = 4, padding=(0.,0.,legend_row_gap,0.))

    Legend(wait_subplot[1, :, Bottom()], symbol_mut, labels_mut, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), colgap = 4, rowgap = 4, padding=(0.,0.,0f0, evo_config.fontsize+legend_row_gap))

    ######################

    ax1 = Axis(pheno_subplot[1,1:4],ygridvisible = false,xgridvisible = false)
    ax2 = Axis(pheno_subplot[1,5],ygridvisible = false,xgridvisible = false)

    m1_mc = reduce(vcat,[tr_mc[1] for (tr_mc,tr) in zip(all_pheno_mut_chf,trajectories)]);
    m1s0_mc = reduce(vcat,reduce(vcat,[tr_mc[2:tr.H0-2] for (tr_mc,tr) in zip(all_pheno_mut_chf,trajectories)]));
    ms0_mc = reduce(vcat,[tr_mc[tr.H0-1] for (tr_mc,tr) in zip(all_pheno_mut_chf,trajectories)]);
    mend_mc = reduce(vcat,reduce(vcat,[tr_mc[tr.H0:end] for (tr_mc,tr) in zip(all_pheno_mut_chf,trajectories)]));

    hm_mc = reduce(hcat,[[count(x->x==i,v)/length(v) for i in [:clb,:crb,:cbb, :mlb,:mrb,:mbb,:neutral, :other]] for v in [m1_mc,m1s0_mc,ms0_mc,mend_mc]])

    all_mchf = reduce(vcat,reduce(vcat,all_pheno_mut_chf))
    all_wait_times = reduce(vcat,map(tr->tr.wait_times[2:end],trajectories));

    av_wait_time = reshape([mean(all_wait_times[all_mchf .== i]) for i in [:clb,:crb,:cbb, :mlb,:mrb,:mbb,:neutral, :other]],(1,8))

    CairoMakie.heatmap!(ax1,1:4,1:8, hm_mc |> transpose)

    CairoMakie.text!(ax1,[Point2f(x, y) for x in 1:4 for y in 1:8],text=
        string.(round.(100 .* vec(hm_mc),digits = 1)) .* "%",
        align = (:center, :center),
        color = ifelse.(hm_mc .< 0.1, :white, :black),fontsize = 10.)

    CairoMakie.heatmap!(ax2,av_wait_time,colormap = :thermal)

    CairoMakie.text!(ax2,[Point2f(1, y) for y in 1:8],text=
        string.(Int.(floor.(vec(av_wait_time)))),
        align = (:center, :center),
        color = ifelse.(av_wait_time .< 5000, :white, :black),fontsize = 10.)

    ax1.xticks = (1:4,[L"n=1",L"1<n<S_0",L"n=S_0",L"n>S_0"])
    ax1.yticks = (1:8,[L"c(L_B)",L"c(R_B)",L"c(LR_B)",L"m(L_B)",L"m(R_B)",L"m(LR_B)",L"\text{Neutral}",L"\text{Other}"])

    hideydecorations!(ax2)

    ax2.xticks = (1:1, [L"\text{Av. wait time}"])

    colgap!(pheno_subplot, Relative(0.01))

    rowgap!(fig.layout, Relative(0.1))
    colgap!(fig.layout, Relative(0.01))

end

function create_weight_edit_summary!(fig,n,trajectories,mutation_operator::MutationOperatorDual,sorted_uep, vertex_top_map, evo_config,color_scheme)

    weight_names_latex = reshape([L"W_{aa}",L"W_{ab}",L"W_{ac}",L"W_{ba}",L"W_{bb}",L"W_{bc}",L"W_{ca}",L"W_{cb}",L"W_{cc}",L"W_{ma}",L"W_{mb}",L"W_{mc}"],(3,4));

    grid_values = Tuple.(findall(ones(3,4) .> 0))

    colors = reverse(palette(:tab10)[1:4])

    ax_wait_list = []

    plot_geno = fig[grid_values[11]...] = GridLayout()

    ax_geno = Axis(plot_geno[1,1],backgroundcolor = (color_scheme[n],1.),title =L"\text{Minimal Stripe Topology}",aspect = DataAspect())

    draw_grn!(ax_geno,vertex_top_map[sorted_uep[n]],evo_config.draw_config,evo_config.node_colors,evo_config.fontsize,true,false)

    norm_type = :pdf

    ###############################

    plot_mut_hist = fig[grid_values[12]...] = GridLayout()

    bins = 50

    min_t_list = [[],[]]
    max_t_list = [[],[]]

    weight_grid_layouts = []

    additive_ax = []
    mult_ax = []

    for type in [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]

        if type[2] == :additive

            min_t = -quantile(mutation_operator.additive_noise_distribution,0.99)
            max_t = quantile(mutation_operator.additive_noise_distribution,0.99)

            mut_noise_dist = mutation_operator.additive_noise_distribution;

            norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];

            if type[1] 
                if n==1
                    ax1 = Axis(plot_mut_hist[1,1],title = L"t<H_{0}")
                    ax2 = Axis(plot_mut_hist[1,2],title = L"t=H_{0}")
                    ax3 = Axis(plot_mut_hist[1,3],title = L"t>H_{0}")
                else
                    ax1 = Axis(plot_mut_hist[1,1])
                    ax2 = Axis(plot_mut_hist[1,2])
                    ax3 = Axis(plot_mut_hist[1,3])
                end

                plot_color = palette(:viridis, 4)[1]

                hidedecorations!(ax1)
            else

                ax1 = Axis(plot_mut_hist[2,1])
                ax2 = Axis(plot_mut_hist[2,2])
                ax3 = Axis(plot_mut_hist[2,3])

                plot_color = palette(:viridis, 4)[2]
            end

            mut_size_r = map(tr->get_mut_size_by_type(tr,type,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

            mut_size = reduce(vcat,mut_size_r)
            
            if length(mut_size) != 0
                CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = plot_color) # new additive
            else
                mut_size = [0.]
            end

            CairoMakie.lines!(ax1,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
            CairoMakie.vlines!(ax1,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth) # new all

            hidedecorations!(ax1)

            mut_size_r = map(tr->get_mut_size_by_type(tr,type,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

            mut_size = reduce(vcat,mut_size_r)

            if length(mut_size) != 0
                CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = plot_color)  # new additive
            else
                mut_size = [0.]
            end

            CairoMakie.lines!(ax2,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
            CairoMakie.vlines!(ax2,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth) # new all

            hidedecorations!(ax2)

            mut_size_r = map(tr->get_mut_size_by_type(tr,type,tr.H0,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

            mut_size = reduce(vcat,mut_size_r)

            if length(mut_size) != 0
                CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = plot_color) # new additive
            else
                mut_size = [0.]
            end

            CairoMakie.lines!(ax3,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
            CairoMakie.vlines!(ax3,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth)

            hidedecorations!(ax3)

            push!(additive_ax,ax1)
            push!(additive_ax,ax2)
            push!(additive_ax,ax3)

        else
            min_t = 0.
            max_t = quantile(mutation_operator.mult_noise_distribution,0.99)

            mut_noise_dist = mutation_operator.mult_noise_distribution;

            norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];

            if type[1] 

                ax1 = Axis(plot_mut_hist[3,1])
                ax2 = Axis(plot_mut_hist[3,2])
                ax3 = Axis(plot_mut_hist[3,3])

                plot_color = palette(:viridis, 4)[3]
                hidedecorations!(ax1)
            else
                ax1 = Axis(plot_mut_hist[4,1])
                ax2 = Axis(plot_mut_hist[4,2])
                ax3 = Axis(plot_mut_hist[4,3])

                plot_color = palette(:viridis, 4)[4]
            end

            mut_size_r = map(tr->get_mut_size_by_type(tr,type,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

            mut_size = abs.(reduce(vcat,mut_size_r))
            
            if length(mut_size) != 0
                CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = plot_color) # new additive
            else
                mut_size = [0.]
            end

            CairoMakie.lines!(ax1,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
            CairoMakie.vlines!(ax1,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth) # new all

            hidedecorations!(ax1)

            mut_size_r = map(tr->get_mut_size_by_type(tr,type,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

            mut_size = abs.(reduce(vcat,mut_size_r))

            if length(mut_size) != 0
                CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = plot_color)  # new additive
            else
                mut_size = [0.]
            end

            CairoMakie.lines!(ax2,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
            CairoMakie.vlines!(ax2,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth) # new all

            hidedecorations!(ax2)

            mut_size_r = map(tr->get_mut_size_by_type(tr,type,tr.H0,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

            mut_size = abs.(reduce(vcat,mut_size_r))

            if length(mut_size) != 0
                CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = plot_color) # new additive
            else
                mut_size = [0.]
            end

            CairoMakie.lines!(ax3,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
            CairoMakie.vlines!(ax3,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth)

            hidedecorations!(ax3)

            push!(mult_ax,ax1)
            push!(mult_ax,ax2)
            push!(mult_ax,ax3)
        end

    end

    linkxaxes!(additive_ax...)
    linkxaxes!(mult_ax...)

    colgap!(plot_mut_hist, 10)
    rowgap!(plot_mut_hist, 10)

    push!(weight_grid_layouts,plot_mut_hist)

    ###############################\ get_mut_size_by_type_and_weight(tr,type,weight_id,1,tr.H0-2)

    bins = 15

    for weight_id in 1:10

        plot_mut_hist = fig[grid_values[weight_id]...] = GridLayout()

        additive_ax = []
        mult_ax = []
    
        for type in [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]

            if type[2] == :additive
    
                min_t = -quantile(mutation_operator.additive_noise_distribution,0.99)
                max_t = quantile(mutation_operator.additive_noise_distribution,0.99)
    
                mut_noise_dist = mutation_operator.additive_noise_distribution;

                norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];
    
                if type[1] 
                    if n==1
                        ax1 = Axis(plot_mut_hist[1,1],title = L"t<H_{0}")
                        ax2 = Axis(plot_mut_hist[1,2],title = L"t=H_{0}")
                        ax3 = Axis(plot_mut_hist[1,3],title = L"t>H_{0}")
                    else
                        ax1 = Axis(plot_mut_hist[1,1])
                        ax2 = Axis(plot_mut_hist[1,2])
                        ax3 = Axis(plot_mut_hist[1,3])
                    end
    
                    plot_color = palette(:viridis, 4)[1]
    
                    hidedecorations!(ax1)
                else
    
                    ax1 = Axis(plot_mut_hist[2,1])
                    ax2 = Axis(plot_mut_hist[2,2])
                    ax3 = Axis(plot_mut_hist[2,3])
    
                    plot_color = palette(:viridis, 4)[2]
                end
    
                mut_size_r = map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))
    
                mut_size = reduce(vcat,mut_size_r)
                
                if length(mut_size) != 0
                    CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = plot_color) # new additive
                else
                    mut_size = [0.]
                end
    
                CairoMakie.lines!(ax1,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
                CairoMakie.vlines!(ax1,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth) # new all
    
                hidedecorations!(ax1)
    
                mut_size_r = map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))
    
                mut_size = reduce(vcat,mut_size_r)
    
                if length(mut_size) != 0
                    CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = plot_color)  # new additive
                else
                    mut_size = [0.]
                end
    
                CairoMakie.lines!(ax2,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
                CairoMakie.vlines!(ax2,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth) # new all
    
                hidedecorations!(ax2)
    
                mut_size_r = map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,tr.H0,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))
    
                mut_size = reduce(vcat,mut_size_r)
    
                if length(mut_size) != 0
                    CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = plot_color) # new additive
                else
                    mut_size = [0.]
                end
    
                CairoMakie.lines!(ax3,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
                CairoMakie.vlines!(ax3,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth)
    
                hidedecorations!(ax3)
    
                push!(additive_ax,ax1)
                push!(additive_ax,ax2)
                push!(additive_ax,ax3)
    
            else
                min_t = 0.
                max_t = quantile(mutation_operator.mult_noise_distribution,0.99)
    
                mut_noise_dist = mutation_operator.mult_noise_distribution;
    
                norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];
    
                if type[1] 

                    ax1 = Axis(plot_mut_hist[3,1])
                    ax2 = Axis(plot_mut_hist[3,2])
                    ax3 = Axis(plot_mut_hist[3,3])
    
                    plot_color = palette(:viridis, 4)[3]
                    hidedecorations!(ax1)
                else
                    ax1 = Axis(plot_mut_hist[4,1])
                    ax2 = Axis(plot_mut_hist[4,2])
                    ax3 = Axis(plot_mut_hist[4,3])

                    plot_color = palette(:viridis, 4)[4]
                end
    
                mut_size_r = map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))
    
                mut_size = abs.(reduce(vcat,mut_size_r))
                
                if length(mut_size) != 0
                    CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = plot_color) # new additive
                else
                    mut_size = [0.]
                end
    
                CairoMakie.lines!(ax1,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
                CairoMakie.vlines!(ax1,1,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth) # new all
    
                hidedecorations!(ax1)
    
                mut_size_r = map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))
    
                mut_size = abs.(reduce(vcat,mut_size_r))
    
                if length(mut_size) != 0
                    CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = plot_color)  # new additive
                else
                    mut_size = [0.]
                end
    
                CairoMakie.lines!(ax2,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
                CairoMakie.vlines!(ax2,1,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth) # new all
    
                hidedecorations!(ax2)
    
                mut_size_r = map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,tr.H0,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))
    
                mut_size = abs.(reduce(vcat,mut_size_r))
    
                if length(mut_size) != 0
                    CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = plot_color) # new additive
                else
                    mut_size = [0.]
                end
    
                CairoMakie.lines!(ax3,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
                CairoMakie.vlines!(ax3,1,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth)
    
                hidedecorations!(ax3)
    
                push!(mult_ax,ax1)
                push!(mult_ax,ax2)
                push!(mult_ax,ax3)
            end
    
        end
    
        linkxaxes!(additive_ax...)
        linkxaxes!(mult_ax...)
    
        colgap!(plot_mut_hist, 10)
        rowgap!(plot_mut_hist, 10)
    
        push!(weight_grid_layouts,plot_mut_hist)
    end

    for (label, layout) in zip(weight_names_latex, weight_grid_layouts[2:end])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = evo_config.fontsize,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right)
    end

    Label(weight_grid_layouts[1][1, 1, TopLeft()], L"\text{All}",
    fontsize = evo_config.fontsize,
    font = :bold,
    padding = (0, 5, 5, 0),
    halign = :right)

    colors = palette(:viridis, 4)

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

#######################

function create_example_traj_fig!(fig_sub,trajectories,sorted_uep,example_mst,tr_choice,ds_config,predict_colors,top_n)

    tr_data = filter(tr->tr.inc_metagraph_vertices[tr.H0] ∈ sorted_uep[example_mst],trajectories);

    example_traj = tr_data[tr_choice];

    ft = example_traj.fitness_traj_tuple[1:example_traj.H0]
    tr_traj = example_traj.topologies[1:example_traj.H0]
    tr_networks = example_traj.geno_traj[1:example_traj.H0]

    tr_we = example_traj.weight_edits[1:example_traj.H0]    
    tr_te = example_traj.top_edits[1:example_traj.H0]    

    tr_prob = vcat(example_traj.gt_label_probabilities[1:example_traj.H0-1,:],reshape([i == example_mst ? 1. : 0. for i in 1:top_n+1],(1,top_n+1)))
    tr_ent = vcat(example_traj.gt_label_entropies[1:example_traj.H0-1],[0.])

    n_networks = length(tr_networks)

    ax_prob_list = []

    for (i,n) in enumerate(tr_networks)

        annotations = string.(round.(reshape(n,(3,4)),digits = 2))

        if i == n_networks
            ax_geno = Axis(fig_sub[4:6,i], backgroundcolor = (predict_colors[example_mst],ds_config_12.color_fade),aspect = DataAspect())
        else
            ax_geno = Axis(fig_sub[4:6,i], backgroundcolor = RGBf(0.98, 0.98, 0.98),aspect = DataAspect())
        end

        draw_grn!(ax_geno,n,ds_config.draw_config,ds_config.node_colors,ds_config.fontsize,annotations,false,false)

    end

    ax_prob = Axis(fig_sub[1:3,1:n_networks], ylabel = L"\text{Probability}",ygridvisible = false,xgridvisible = false)
    ax_pro = Axis(fig_sub[1:3,1:n_networks],yaxisposition = :right, yticklabelcolor = :green)

    hidespines!(ax_pro)
    hideydecorations!(ax_pro,label = false,ticklabels = false,ticks = false,minorticks = false)
    hidexdecorations!(ax_pro)

    px = reduce(vcat,[[n for i in 1:size(tr_prob,2)] for n in 1:size(tr_prob,1)])
    pgrp = reduce(vcat,[[i for i in 1:size(tr_prob,2)] for n in 1:size(tr_prob,1)])

    CairoMakie.barplot!(ax_prob,px,vec(tr_prob |> transpose |> collect),color = [predict_colors[i] for i in pgrp], dodge = pgrp, gap = 0.55)

    ent_l = CairoMakie.lines!(ax_prob,tr_ent)
    ent_s = CairoMakie.scatter!(ax_prob,tr_ent, marker = 'x')

    we_l = CairoMakie.lines!(ax_pro,tr_we, color = :green)
    we_s = CairoMakie.scatter!(ax_pro,tr_we, marker = 'x', color = :green)

    te_l = CairoMakie.lines!(ax_pro,tr_te, color = :green, linestyle = "--")
    te_s = CairoMakie.scatter!(ax_pro,tr_te, marker = 'x', color = :green)

    ax_prob.xticks = (1:5,[L"p^{(i)}(1)",L"p^{(i)}(S_0-3)",L"p^{(i)}(S_0-2)",L"p^{(i)}(S_0-1)",L"p^{(i)}(S_0)"])

    CairoMakie.ylims!(ax_prob,0.,1.)

    ax_pro.yticks = ([i for i in 1:maximum(tr_we) if i%2 == 1],string.([i for i in 1:maximum(tr_we) if i%2 == 1]))

    Legend(fig_sub[7, 1:n_networks],[[ent_l,ent_s],[te_l,te_s],[we_l,we_s]],[L"H(p^{(i)}(n))",L"\text{Cumul. topology changes}",L"\text{Cumul. interaction changes}"],  colgap = 4, padding = (0.,0.,0.,0.),framevisible=false, orientation = :horizontal,nbanks = 1)

    linkxaxes!(ax_pro,ax_prob)

    rowgap!(fig_sub, Relative(0.01))
    colgap!(fig_sub, Relative(0.02))

end

function create_pred_accuracy_fig!(fig_sub,streak,label_streak,predict_colors,top_n)

    ax_list = []
    ax_list_h = []

    vl = []

    for n in 1:top_n+1

        if n == 1
            ax = Axis(fig_sub[1,n],ylabel  = L"\text{C.D.F}",xgridvisible = false,ygridvisible = false)
        elseif n == 3
            ax = Axis(fig_sub[1,n],xlabel = L"\text{Earliest correct prediction, as % of cumulative mutations at } S_0",xgridvisible = false,ygridvisible = false)
            hideydecorations!(ax)
        else
            ax = Axis(fig_sub[1,n],xgridvisible = false,ygridvisible = false)
            hideydecorations!(ax)
        end

        if n == top_n+1
            ax_hist = Axis(fig_sub[1,n],ylabel  = L"\text{P.D.F}",xgridvisible = false,ygridvisible = false,yaxisposition = :right)
            hideydecorations!(ax_hist,label = false,ticklabels = false,ticks = false,minorticks = false)
            hidexdecorations!(ax_hist)
        else
            ax_hist = Axis(fig_sub[1,n],xgridvisible = false,ygridvisible = false)
            hidespines!(ax_hist)
            hideydecorations!(ax_hist)
            hidexdecorations!(ax_hist)
        end

        n_trials = sum(label_streak .== n)

        cum_streak = [sum(streak[label_streak .== n] .<= i)/n_trials for i in 0:0.1:1-eps()];

        CairoMakie.hist!(ax_hist,streak[label_streak .== n], nbins = 100, color = predict_colors[n], normalization = :probability)

        l1 = CairoMakie.lines!(ax,0:0.1:1-eps(),cum_streak,color = predict_colors[n], linestyle = :dash)

        push!(ax_list,ax)
        push!(ax_list_h,ax_hist)

        push!(vl,l1)

        ax.xticks = ([0.,0.5,1.], ["0","0.5","1"])

        CairoMakie.ylims!(ax,0.,0.75)

        # axislegend(ax,[l1],[L"\text{CDF}"], position = :lt)

    end

    linkyaxes!(ax_list...)
    linkyaxes!(ax_list_h...)

    colgap!(fig_sub, Relative(0.02))
    rowgap!(fig_sub, Relative(0.02))
end

function find_candidates(founder,all_target_networks,contin_target,top_n,model_gtl,test_indices,pred_grid,prob_grid,fitness_grid,sample_points,min_prob)

    contin_choices = []
    contin_prob = []

    founder_fitness = fitness_function(founder.phenotype)

    sample_grid = [[w1,w2] for w1 in sample_points, w2 in sample_points]

    sample_grid_v = reduce(hcat,vec(sample_grid));

    for n in 1:top_n+1

        ########################### Grids
    
        resultant_network = all_target_networks[n][2]
        resultant_network_m = reshape(resultant_network,(3,4))
    
        start_network = all_target_networks[n][1]
        start_network_m = reshape(start_network,(3,4))
    
        mut_id = Tuple.(findall(x->x!=0,reshape(all_target_networks[n][2],(3,4)) .!= reshape(all_target_networks[n][1],(3,4))))
    
        t1 = mut_id[1]
        t2 = mut_id[2]
    
        n1 = findall(x->x == t1,test_indices)[1]
        n2 = findall(x->x == t2,test_indices)[1]
    
        ###############
    
        contin_choice = first(Tuple(get_max_prob_point(n1,n2,pred_grid,prob_grid,fitness_grid,contin_target[n],min_prob,founder_fitness)))
    
        w1 = sample_grid_v[1,contin_choice]
        w2 = sample_grid_v[2,contin_choice]
    
        contin_w = set_weight(start_network_m,t1,t2,w1,w2)
        contin_wm = reshape(contin_w,(3,4))
    
        gt_dtrain = xgboost.DMatrix(contin_w[1:10] |> transpose |> collect, feature_names = weight_names)
        gt_label_probabilities = model_gtl.predict(gt_dtrain)
    
        push!(contin_choices,contin_wm)
        push!(contin_prob,gt_label_probabilities)
    end

    return contin_choices,contin_prob
end

function create_fm_result_fig!(fig_result_sub,fig_fl_sub,plot_save_dir,founder,all_target_networks,all_target_prob,vim_trajectories,contin_choices,top_n,mut_q,predict_label_to_vertex,draw_config_12,fs12_default,node_colors,fontsize_pub,predict_colors,save_top_fig = true)

    ax_list = []

    vpreds = []
    vbars = []

    fitness_prob_all = []

    founder_fitness = fitness_function(founder.phenotype)

    for n in 1:top_n+1

        mut_id = Tuple.(findall(x->x!=0,reshape(all_target_networks[n][2],(3,4)) .!= reshape(all_target_networks[n][1],(3,4))))

        t1 = mut_id[1]
        t2 = mut_id[2]

        n1 = findall(x->x == t1,test_indices)[1]
        n2 = findall(x->x == t2,test_indices)[1]

        # prob_fv = [log10(fixation_probability_kim(0.,f,β[1],β[2])) for f in fitness_grid[n1,n2,:] .- founder_fitness[2]]

        prob_fv = log10.([fixation_probability_kim(f[1] - founder_fitness[1],f[2] - founder_fitness[2],β[1],β[2]) for f in fitness_grid[n1,n2,:]])

        push!(fitness_prob_all,prob_fv)
    end

    prob_v = reduce(vcat,fitness_prob_all)
    min_prob = minimum(prob_v[prob_v .> log10(2e-53)])

    for n in 1:top_n+1

        if n == 1
            ax = Axis(fig_result_sub[1,n],ylabel  = L"\text{Freq}", xlabel = L"M^{(i)}_{N_i}")
            hidexdecorations!(ax)
        else
            ax = Axis(fig_result_sub[1,n],xlabel = L"M^{(i)}_{N_i}")
            hideydecorations!(ax)
            hidexdecorations!(ax)
        end

        axm_fit = Axis(fig_fl_sub[1,n])
        axm_pred = Axis(fig_fl_sub[2,n])

        #################

        pop = vim_trajectories[n]

        incl_count_all = []
        
        for k in 1:top_n+1
            incl_count = sum(map(tr->tr.inc_metagraph_vertices == predict_label_to_vertex[k],pop))
            push!(incl_count_all,incl_count)
        end

        p_dist = incl_count_all ./ sum(incl_count_all)

        push!(other,incl_count_all)

        p_dist_x = [i for i in 1:top_n+1]

        dodge_current = [1 for _ in 1:top_n+1]

        dodge_H0 = [2 for _ in 1:top_n+1]

        x = vcat(p_dist_x,p_dist_x)

        color = vcat([(predict_colors[i],1.) for i in p_dist_x],[:white for i in p_dist_x])

        stroke_color = [(predict_colors[i],1.) for i in x]

        dodge = vcat(dodge_current,dodge_H0)
        y = vcat(p_dist,null_H0_dist)

        vbar = CairoMakie.barplot!(ax,x,y,dodge = dodge,color = color, strokecolor = stroke_color, strokewidth = 1., gap = 0.2)

        vpred = CairoMakie.scatter!(ax,[1,2,3,4,5], all_target_prob[n], marker = 'x', color = :cyan, markersize = 10.)

        push!(vbars,vbar)
        push!(vpreds,vpred)

        CairoMakie.hidedecorations!(ax,label = false,ticklabels = false,ticks = false,minorticks = false)

        axis_label = vcat(collect(string.(1:top_n)),["o"])
        
        ax.xticks = (1:top_n+1,axis_label)

        CairoMakie.ylims!(ax,0.,1.)

        push!(ax_list,ax)

        start_network = all_target_networks[n][1]
        resultant_network = all_target_networks[n][2]

        if save_top_fig
            fig_spare_m1 = CairoMakie.Figure(resolution = (100,100))
            ax_geno_m = Axis(fig_spare_m1[1,1],aspect = DataAspect())

            draw_grn_mutant!(ax_geno_m,start_network,resultant_network,draw_config_12,fs12_default,node_colors,fontsize_pub,predict_colors[n],false,false)

            cond_save(plotsdirx(plot_save_dir,"Di-MST" * string(n) * ".pdf"),fig_spare_m1,true)
        end

        ########################### Grids
        
        resultant_network = all_target_networks[n][2]
        resultant_network_m = reshape(resultant_network,(3,4))

        start_network = all_target_networks[n][1]
        start_network_m = reshape(start_network,(3,4))

        mut_id = Tuple.(findall(x->x!=0,reshape(all_target_networks[n][2],(3,4)) .!= reshape(all_target_networks[n][1],(3,4))))

        t1 = mut_id[1]
        t2 = mut_id[2]

        n1 = findall(x->x == t1,test_indices)[1]
        n2 = findall(x->x == t2,test_indices)[1]

        ##############

        ylabelwl = Label(fig_fl_sub[1:2,n,Left()], weight_names_latex_m[t2...],  padding = (0.,2.,0.,0.))
        ylabelwl = Label(fig_fl_sub[1,n,Bottom()], "",  padding = (0.,0.,5.,0.))

        sn1 = start_network_m[t1...]
        sn2 = start_network_m[t2...]

        rn1 = resultant_network_m[t1...]
        rn2 = resultant_network_m[t2...]
        
        rnc1 = contin_choices[n][t1...]
        rnc2 = contin_choices[n][t2...]

        fit_v = log10.([fixation_probability_kim(f[1] - founder_fitness[1],f[2] - founder_fitness[2],β[1],β[2]) for f in fitness_grid[n1,n2,:]])

        hm = CairoMakie.heatmap!(axm_fit,sample_points,sample_points,reshape(fit_v,(N_sample,N_sample)),colormap = cgrad(:viridis), colorrange = (min_prob,0.), lowclip = :black)

        if n == 1 
            cb = Colorbar(fig_fl_sub[1:2, top_n+2], hm;
            minorticksvisible=true,width = 5.,tickformat = custom_log_formatter,ticklabelsize = 9.
                )
        end

        max_ent = maximum(1 .- entropy_grid[n1,n2,:])
        min_ent = minimum(1 .- entropy_grid[n1,n2,:])

        pred_vc = [log10(fixation_probability_kim(f - founder_fitness[2],β[1],β[2])) >  min_prob ? (predict_colors[c[1]],c[2]) : :black for (c,f) in zip(zip(pred_grid[n1,n2,:],prob_grid[n1,n2,:]),fitness_grid[n1,n2,:])]

        CairoMakie.scatter!(axm_pred,sample_grid_v, color = pred_vc, markersize = 2.)

        CairoMakie.scatter!(axm_pred,reshape([sn1,sn2],(2,1)), color = :yellow, markersize = 3.)
        CairoMakie.scatter!(axm_pred,reshape([rn1,rn2],(2,1)), color = :cyan, markersize = 3.)

        arrow_start = [Point2f(sn1,sn2)]
        arrow_dir = [Vec2f(rn1-sn1,rn2-sn2)]

        arrow_dir_c = [Vec2f(rnc1-sn1,rnc2-sn2)]

        CairoMakie.arrows!(axm_pred,arrow_start,arrow_dir,arrowsize = 5., color = :cyan)
        CairoMakie.arrows!(axm_fit,arrow_start,arrow_dir,arrowsize = 5., color = :cyan)

        CairoMakie.arrows!(axm_pred,arrow_start,arrow_dir_c,arrowsize = 5., color = :yellow)
        CairoMakie.arrows!(axm_fit,arrow_start,arrow_dir_c,arrowsize = 5., color = :yellow)

        xmin = ceil(Int,sn1 - quantile(mutation_op.additive_noise_distribution,mut_q))
        xmax = ceil(Int,sn1 + quantile(mutation_op.additive_noise_distribution,mut_q))

        ymin = ceil(Int,sn2 - quantile(mutation_op.additive_noise_distribution,mut_q))
        ymax = ceil(Int,sn2 + quantile(mutation_op.additive_noise_distribution,mut_q))

        axm_pred.xlabel = weight_names_latex_m[t1...]

        CairoMakie.xlims!(axm_pred,xmin,xmax)
        CairoMakie.ylims!(axm_pred,ymin,ymax)

        hidexdecorations!(axm_pred,label = false)
        hideydecorations!(axm_pred,label = false)

        CairoMakie.xlims!(axm_fit,xmin,xmax)
        CairoMakie.ylims!(axm_fit,ymin,ymax)

        hidexdecorations!(axm_fit,label = false)
        hideydecorations!(axm_fit,label = false)

        vlines!(axm_pred,0., linestyle = :dash, color = :black, linewidth = 1.)
        hlines!(axm_pred,0., linestyle = :dash, color = :black, linewidth = 1.)

        vlines!(axm_fit,0., linestyle = :dash, color = :black, linewidth = 1.)
        hlines!(axm_fit,0., linestyle = :dash, color = :black, linewidth = 1.)
    end

    linkyaxes!(ax_list...)

    colgap!(fig_fl_sub, Relative(0.01))
    rowgap!(fig_fl_sub, Relative(0.005))

    colgap!(fig_result_sub, Relative(0.02))
    rowgap!(fig_result_sub, Relative(0.075))
end

function create_fm_result_partial_fig!(fig_result_sub,fig_fl_sub,founder,all_target_networks,top_n,sample_points)

    ax_list = []

    vpreds = []
    vbars = []

    fitness_prob_all = []

    founder_fitness = fitness_function(founder.phenotype)

    sample_grid = [[w1,w2] for w1 in sample_points, w2 in sample_points]

    sample_grid_v = reduce(hcat,vec(sample_grid));

    for n in 1:top_n+1

        mut_id = Tuple.(findall(x->x!=0,reshape(all_target_networks[n][2],(3,4)) .!= reshape(all_target_networks[n][1],(3,4))))

        t1 = mut_id[1]
        t2 = mut_id[2]

        n1 = findall(x->x == t1,test_indices)[1]
        n2 = findall(x->x == t2,test_indices)[1]

        # prob_fv = [log10(fixation_probability_kim(0.,f,β[1],β[2])) for f in fitness_grid[n1,n2,:] .- founder_fitness[2]]

        prob_fv = log10.([fixation_probability_kim(f[1] - founder_fitness[1],f[2] - founder_fitness[2],β[1],β[2]) for f in fitness_grid[n1,n2,:]])

        push!(fitness_prob_all,prob_fv)
    end

    prob_v = reduce(vcat,fitness_prob_all)
    min_prob = minimum(prob_v[prob_v .> log10(2e-53)])

    for n in 1:top_n+1

        if n == 1
            ax = Axis(fig_result_sub[1,n],ylabel  = L"\text{Freq}", xlabel = L"M^{(i)}_{N_i}")
            hidexdecorations!(ax)
        else
            ax = Axis(fig_result_sub[1,n],xlabel = L"M^{(i)}_{N_i}")
            hideydecorations!(ax)
            hidexdecorations!(ax)
        end

        axm_fit = Axis(fig_fl_sub[1,n])
        axm_pred = Axis(fig_fl_sub[2,n])

        #################

        ########################### Grids
        
        resultant_network = all_target_networks[n][2]
        resultant_network_m = reshape(resultant_network,(3,4))

        start_network = all_target_networks[n][1]
        start_network_m = reshape(start_network,(3,4))

        mut_id = Tuple.(findall(x->x!=0,reshape(all_target_networks[n][2],(3,4)) .!= reshape(all_target_networks[n][1],(3,4))))

        t1 = mut_id[1]
        t2 = mut_id[2]

        n1 = findall(x->x == t1,test_indices)[1]
        n2 = findall(x->x == t2,test_indices)[1]

        ##############

        ylabelwl = Label(fig_fl_sub[1:2,n,Left()], weight_names_latex_m[t2...],  padding = (0.,2.,0.,0.))
        ylabelwl = Label(fig_fl_sub[1,n,Bottom()], "",  padding = (0.,0.,5.,0.))

        sn1 = start_network_m[t1...]
        sn2 = start_network_m[t2...]

        rn1 = resultant_network_m[t1...]
        rn2 = resultant_network_m[t2...]

        fit_v = log10.([fixation_probability_kim(f[1] - founder_fitness[1],f[2] - founder_fitness[2],β[1],β[2]) for f in fitness_grid[n1,n2,:]])

        hm = CairoMakie.heatmap!(axm_fit,sample_points,sample_points,reshape(fit_v,(N_sample,N_sample)),colormap = cgrad(:viridis), colorrange = (min_prob,0.), lowclip = :black)

        if n == 1 
            cb = Colorbar(fig_fl_sub[1:2, top_n+2], hm;
            minorticksvisible=true,width = 5.,tickformat = custom_log_formatter,ticklabelsize = 9.
                )
        end

        max_ent = maximum(1 .- entropy_grid[n1,n2,:])
        min_ent = minimum(1 .- entropy_grid[n1,n2,:])

        pred_vc = [log10(fixation_probability_kim(f[1] - founder_fitness[1],f[2] - founder_fitness[2],β[1],β[2])) >  min_prob ? (predict_colors[c[1]],c[2]) : :black for (c,f) in zip(zip(pred_grid[n1,n2,:],prob_grid[n1,n2,:]),fitness_grid[n1,n2,:])]

        CairoMakie.scatter!(axm_pred,sample_grid_v, color = pred_vc, markersize = 2.)

        CairoMakie.scatter!(axm_pred,reshape([sn1,sn2],(2,1)), color = :yellow, markersize = 3.)
        CairoMakie.scatter!(axm_pred,reshape([rn1,rn2],(2,1)), color = :cyan, markersize = 3.)

        arrow_start = [Point2f(sn1,sn2)]
        arrow_dir = [Vec2f(rn1-sn1,rn2-sn2)]

        CairoMakie.arrows!(axm_pred,arrow_start,arrow_dir,arrowsize = 5., color = :cyan)
        CairoMakie.arrows!(axm_fit,arrow_start,arrow_dir,arrowsize = 5., color = :cyan)

        xmin = ceil(Int,sn1 - quantile(mutation_op.additive_noise_distribution,mut_q))
        xmax = ceil(Int,sn1 + quantile(mutation_op.additive_noise_distribution,mut_q))

        ymin = ceil(Int,sn2 - quantile(mutation_op.additive_noise_distribution,mut_q))
        ymax = ceil(Int,sn2 + quantile(mutation_op.additive_noise_distribution,mut_q))

        axm_pred.xlabel = weight_names_latex_m[t1...]

        CairoMakie.xlims!(axm_pred,xmin,xmax)
        CairoMakie.ylims!(axm_pred,ymin,ymax)

        hidexdecorations!(axm_pred,label = false)
        hideydecorations!(axm_pred,label = false)

        CairoMakie.xlims!(axm_fit,xmin,xmax)
        CairoMakie.ylims!(axm_fit,ymin,ymax)

        hidexdecorations!(axm_fit,label = false)
        hideydecorations!(axm_fit,label = false)

        vlines!(axm_pred,0., linestyle = :dash, color = :black, linewidth = 1.)
        hlines!(axm_pred,0., linestyle = :dash, color = :black, linewidth = 1.)

        vlines!(axm_fit,0., linestyle = :dash, color = :black, linewidth = 1.)
        hlines!(axm_fit,0., linestyle = :dash, color = :black, linewidth = 1.)
    end

    # linkyaxes!(ax_list...)

    colgap!(fig_fl_sub, Relative(0.01))
    rowgap!(fig_fl_sub, Relative(0.005))

    colgap!(fig_result_sub, Relative(0.02))
    rowgap!(fig_result_sub, Relative(0.075))
end

function create_fm_result_fig!(fig_result_sub,plot_save_dir,vim_trajectories,contin_choices,top_n,mut_q,predict_label_to_vertex,draw_config_12,fs12_default,node_colors,fontsize_pub,predict_colors,save_top_fig = true)

    ax_list = []

    for n in 1:top_n+1

        fig_spare_m1 = CairoMakie.Figure(resolution = (100,100))
        ax_geno_m = Axis(fig_spare_m1[1,1],aspect = DataAspect())

        if n == 1
            ax = Axis(fig_result_sub[2,n],ylabel  = L"\text{Freq}", xlabel = L"M^{(i)}_{N_i}")
            hidexdecorations!(ax)
        else
            ax = Axis(fig_result_sub[2,n],xlabel = L"M^{(i)}_{N_i}")
            hideydecorations!(ax)
            hidexdecorations!(ax)
        end


        if n == top_n+1
            p_dist = [0.,0.,0.,0.,0.]
        else
            pop = vim_trajectories[n]

            incl_count_all = []
        
            for k in 1:top_n+1
                incl_count = sum(map(tr->tr.inc_metagraph_vertices == predict_label_to_vertex[k],pop))
                push!(incl_count_all,incl_count)
            end

            p_dist = incl_count_all ./ sum(incl_count_all)
        end

        p_dist_x = [i for i in 1:top_n+1]

        dodge_current = [1 for _ in 1:top_n+1]

        dodge_H0 = [2 for _ in 1:top_n+1]

        x = vcat(p_dist_x,p_dist_x)

        color = vcat([(predict_colors[i],1.) for i in p_dist_x],[:white for i in p_dist_x])

        stroke_color = [(predict_colors[i],1.) for i in x]

        dodge = vcat(dodge_current,dodge_H0)
        y = vcat(p_dist,null_H0_dist)

        vbar = CairoMakie.barplot!(ax,x,y,dodge = dodge,color = color, strokecolor = stroke_color, strokewidth = 1., gap = 0.2)

        print(n)

        vpred = CairoMakie.scatter!(ax,[1,2,3,4,5], contin_prob[n][1,:], marker = 'x', color = :yellow, markersize = 10.)

        CairoMakie.hidedecorations!(ax,label = false,ticklabels = false,ticks = false,minorticks = false)

        # print(n)

        axis_label = vcat(collect(string.(1:top_n)),["o"])
        
        ax.xticks = (1:top_n+1,axis_label)

        push!(ax_list,ax)

        start_network = all_target_networks[n][1]
        resultant_network = vec(contin_choices[n])

        draw_grn_mutant!(ax_geno_m,start_network,resultant_network,draw_config_12,fs12_default,node_colors,fontsize_pub,predict_colors[contin_target[n]],false,false)

        cond_save(plotsdirx(plot_save_dir,"Dii-MST" * string(n) * ".pdf"),fig_spare_m1,true)
    end

    colgap!(fig_result_sub, Relative(0.02))
    rowgap!(fig_result_sub, Relative(0.075))

    linkyaxes!(ax_list...)
end

function create_ent_hist_fig!(fig_hist_sub,fig_pie_sub,trajectories_p,vertex_to_predict_label,top_n,predict_colors,evo_config_12)

    ax = Axis(fig_hist_sub[1,1],ygridvisible = false,xgridvisible = false, xlabel = L"H(p^{(i)}(1))", ylabel = L"\text{Freq}")

    traj_epi_pop = filter(tr->(length(tr.mutant_info[1].weight_id) == 2) & (tr.gt_label_predictions[2] .!= -1),trajectories_p)

    all_first_entropies = map(tr->tr.gt_label_entropies[2],traj_epi_pop)

    all_first_pred = [tr.gt_label_predictions[2] == -1 ? 5 : vertex_to_predict_label[tr.gt_label_predictions[2]] for tr in traj_epi_pop]

    first_epi = map(tr->tr.epistasis[2],traj_epi_pop)

    CairoMakie.hist!(ax,all_first_entropies,bins = 50,color = (:grey,0.5), normalization = :probability)

    bin_boundaries = [0.6,0.65,0.7]
    bin_cp = [0.3,0.625,0.675,0.85]

    CairoMakie.vlines!(ax,bin_boundaries,linestyle = :dash, color = :black)

    h_first_ent = fit(Histogram, all_first_entropies, bin_boundaries; closed = :left) 

    first_ent_bins = map(f->StatsBase.binindex(h_first_ent, f),all_first_entropies);

    for i in 0:length(bin_boundaries)

        ax_pred = Axis(fig_pie_sub[1,i+1]) 

        pred_prop = calculate_pred_proportion(all_first_pred[first_ent_bins .== i],top_n+1)

        CairoMakie.pie!(ax_pred,pred_prop,radius = evo_config_12.pie_radius,color = predict_colors,
        inner_radius = evo_config_12.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config_12.pie_strokewidth)

        CairoMakie.text!(ax_pred,0.,0.,text = [string(i+1)], fontsize = 10., align = (:center,:center))

        CairoMakie.hidedecorations!(ax_pred)

    end

    Label(fig_pie_sub[1,:,Bottom()], L"\text{argmax} p^{(i)}(1)", padding = (0.,0.,40.,0.))

    CairoMakie.text!(ax,bin_cp,[0.05 for _ in bin_cp],text = [string(i) for i in 1:length(bin_cp)], fontsize = 10.,align = (:center,:center))

    rowgap!(fig_hist_sub, Relative(0.02))
    colgap!(fig_hist_sub, Relative(0.02))

    rowgap!(fig_pie_sub, Relative(0.02))
    colgap!(fig_pie_sub, Relative(0.02))
end

function create_full_prediction_summary!(fig,trajectories_p,sorted_uep,top_n,contingency_data,contin_choices,predict_label_to_vertex,vertex_to_predict_label,mut_q,predict_colors,fs_default,evo_config,ds_config,fontsize)

    fig_ex = fig[1:4,1:top_n+2] = GridLayout()

    ######### Example

    example_mst = 4

    tr_choice = 3404

    create_example_traj_fig!(fig_ex,trajectories_p,sorted_uep,example_mst,tr_choice,ds_config,predict_colors,top_n)

    ############# Pred accuracy

    fig_n  = fig[5,1:top_n+2] = GridLayout()

    create_pred_accuracy_fig!(fig_n,streak,label_streak,predict_colors,top_n)

    ############ Contin Double

    fig_fit = fig[6:7,1:top_n+2] = GridLayout()
    fig_d = fig[8:10,1:top_n+2] = GridLayout()

    create_fm_result_fig!(fig_d,fig_fit,plot_save_dir,founder,all_target_networks,all_target_prob,vim_trajectories,contin_choices,top_n,mut_q,predict_label_to_vertex,ds_config.draw_config,fs_default,ds_config.node_colors,fontsize,predict_colors,true)

    vim_trajectories = contingency_data["Single"][1]

    create_fm_result_fig!(fig_d,plot_save_dir,vim_trajectories,contin_choices,top_n,mut_q,predict_label_to_vertex,ds_config.draw_config,fs_default,ds_config.node_colors,fontsize,predict_colors,true)

    # ############ works

    fig_hist = fig[11,1:3] = GridLayout()
    fig_pie = fig[11,4:6] = GridLayout()

    create_ent_hist_fig!(fig_hist,fig_pie,trajectories_p,vertex_to_predict_label,top_n,predict_colors,evo_config)

    ############

    fig_l = fig[12,1:top_n+2] = GridLayout()

    Legend(fig_l[1, 1],[vbars[end],vpreds[1],vl[1]],[L"\text{Solid bars = realised distributions}", L"\text{Predicted outcome probabilities}",L"\text{Empirical CDF}"],framevisible=false, nbanks = 2)

    for (label, layout) in zip(["A", "B", "C","D", "E"], [fig_ex,fig_n,fig_fit,fig_d,fig_hist])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = ds_config.caption_fontsize,
            font = :bold,
            padding = (0,ds_config.caption_padding, ds_config.caption_padding, 0),
            halign = :right)
    end

    colgap!(fig.layout, Relative(0.025))
    rowgap!(fig.layout, Relative(0.005))

end


function create_full_prediction_summary_nr!(fig,trajectories_p,sorted_uep,top_n,predict_colors,sample_points,ds_config)

    fig_ex = fig[1:4,1:top_n+2] = GridLayout()

    ######### Example

    example_mst = 4

    tr_choice = 3404

    create_example_traj_fig!(fig_ex,trajectories_p,sorted_uep,example_mst,tr_choice,ds_config,predict_colors,top_n)

    ############# Pred accuracy

    fig_n  = fig[5,1:top_n+2] = GridLayout()

    create_pred_accuracy_fig!(fig_n,streak,label_streak,predict_colors,top_n)

    ############ Contin Double

    fig_fit = fig[6:7,1:top_n+2] = GridLayout()
    fig_d = fig[8:10,1:top_n+2] = GridLayout()

    create_fm_result_partial_fig!(fig_d,fig_fit,founder,all_target_networks,top_n,sample_points)

    # ############ works

    fig_hist = fig[11,1:3] = GridLayout()
    fig_pie = fig[11,4:6] = GridLayout()

    ############

    fig_l = fig[12,1:top_n+2] = GridLayout()

    colgap!(fig.layout, Relative(0.025))
    rowgap!(fig.layout, Relative(0.005))

end


