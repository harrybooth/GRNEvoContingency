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

function plot_dynamical_summary_parents!(fig,trajectories,embedding,top_n,minimal_motif_count,sorted_uep,sorted_counts_uep,mst_conf_int,end_parents,vertex_top_map,example_mst,tr_choice,ds_config)

    trajectories_p_d = filter(tr->tr.inc_metagraph_parents[end] ∈ sorted_uep[1:top_n],trajectories);

    mo_umap = fig[1:4, 1:4] = GridLayout()

    ex1 = fig[1:3, 5:8] = GridLayout()
    rmh0 = fig[4, 5:8] = GridLayout()

    top_n_dict = Dict(v_id=>pos for (pos,v_id) in enumerate(sorted_uep[1:top_n]))

    #### Motif Distribution

    color_sorted_counts_uep = [i <= top_n ? ds_config.color_scheme[i] : :grey for i in 1:length(sorted_counts_uep)]

    view_sorted_uep_id = sorted_counts_uep .> minimal_motif_count

    # other = mean(sorted_counts_uep[.!view_sorted_uep_id])

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

    tr_data = filter(tr->tr.inc_metagraph_parents[tr.H0] ∈ sorted_uep[example_mst],trajectories);
    tr_data_id = findall(tr->tr.inc_metagraph_parents[tr.H0] ∈ sorted_uep[example_mst],trajectories)

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

        pop = filter(tr->tr.inc_metagraph_parents[end] == sorted_uep[n],trajectories_p_d)

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

        epi_counts = reduce(vcat,map(tr->tr.epistasis[1:tr.H0-2],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

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

        epi_counts = map(tr->tr.epistasis[tr.H0-1],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

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

        epi_counts = reduce(vcat,map(tr->tr.epistasis[tr.H0:end],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

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

        epi_counts = reduce(vcat,map(tr->tr.epistasis[1:tr.H0-2],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

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

        epi_counts = map(tr->tr.epistasis[tr.H0-1],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

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

        epi_counts = reduce(vcat,map(tr->tr.epistasis[tr.H0:end],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

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

function create_evo_summary_portrait!(fig,trajectories,top_n,mutation_operator::Union{MutationOperatorDual,MutationOperatorUniform},sorted_uep, vertex_top_map,wait_time_summary,evo_config)

    all_wait_times = reduce(hcat,[average_wait_time(tr) for tr in trajectories]);

    ax_wait_list = []

    ax_wait_list = []
    ax_wait_2_list = []

    wt_l_list = []
    wt_s_list = []

    min_t_u = -mutation_operator.max_w
    max_t_u = mutation_operator.max_w

    for n in 1:top_n

        plot_geno = fig[n, 1] = GridLayout()
        plot_wait = fig[n, 2:4] = GridLayout()
        plot_mut_hist = fig[n, 5:9] = GridLayout()

        if n==1
            ax_geno = Axis(plot_geno[1,1],backgroundcolor = (evo_config.color_scheme[n],evo_config.color_fade),title =L"M^{(i)}_{N_i}",aspect = DataAspect())
        else
            ax_geno = Axis(plot_geno[1,1],backgroundcolor = (evo_config.color_scheme[n],evo_config.color_fade),aspect = DataAspect())
        end

        ax_geno_s = Axis(plot_geno[2,1],backgroundcolor = :transparent,
        leftspinevisible = false,
        rightspinevisible = false,
        bottomspinevisible = false,
        topspinevisible = false,
        xticklabelsvisible = false, 
        yticklabelsvisible = false,
        xgridcolor = :transparent,
        ygridcolor = :transparent,
        xminorticksvisible = false,
        yminorticksvisible = false,
        xticksvisible = false,
        yticksvisible = false,
        xautolimitmargin = (0.0,0.0),
        yautolimitmargin = (0.0,0.0))

        hidedecorations!(ax_geno_s)

        draw_grn!(ax_geno,vertex_top_map[sorted_uep[n]],evo_config.draw_config,evo_config.node_colors,evo_config.fontsize,false,false)

        ax_wait = Axis(plot_wait[2,1:3],yticklabelsize = 0.8*evo_config.fontsize,yaxisposition = :right, xticklabelsvisible = false, yticksize= 0.25*evo_config.fontsize)

        ax_wait_2 = Axis(plot_wait[2,1:3], yticklabelcolor = :red,yscale = log10,yticklabelsize = 0.8*evo_config.fontsize,yticksize= 0.25*evo_config.fontsize)

        hidespines!(ax_wait_2 )
        hideydecorations!(ax_wait_2,label = false,ticklabels = false,ticks = false,minorticks = false)
        hidexdecorations!(ax_wait_2)

        # #############################

        mut_type_prop_all = []
        mut_type_time_labels = []
        mut_type_labels = []

        # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,1,tr.H0-2), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n]) & (tr.H0-2 > 0),trajectories_p))

        mut_type_prop = map(tr->calculate_mut_type_count(get_mut_type(tr,1,tr.H0-2), [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n]) & (tr.H0-2 > 0),trajectories))

        mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

        push!(mut_type_prop_all,mut_type_prop_av)
        push!(mut_type_labels, [1,2,3,4])
        push!(mut_type_time_labels,[1,1,1,1])

        # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0-1,tr.H0-1), [:existing,:new,:del]),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n] ,trajectories_p))
        mut_type_prop = map(tr->calculate_mut_type_count(get_mut_type(tr,tr.H0-1,tr.H0-1), [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n] ,trajectories))

        mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

        push!(mut_type_prop_all,mut_type_prop_av)
        push!(mut_type_labels, [1,2,3,4])
        push!(mut_type_time_labels,[2,2,2,2])

        # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0,length(tr.topologies)-1), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n])  & (tr.H0 < length(tr.topologies)),trajectories_p))

        mut_type_prop = map(tr->calculate_mut_type_count(get_mut_type(tr,tr.H0,length(tr.topologies)-1), [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n])  & (tr.H0 < length(tr.topologies)),trajectories))

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

        ax_epi_lH0 = Axis(plot_wait[1,1],title = L"t<H_{0}")

        epi_counts = reduce(vcat,map(tr->tr.epistasis[1:tr.H0-2],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_lH0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_lH0)

        ax_epi_H0 = Axis(plot_wait[1,2],title = L"t=H_{0}")
 
        epi_counts = map(tr->tr.epistasis[tr.H0-1],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_H0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_H0)

        ax_epi_uH0 = Axis(plot_wait[1,3],title = L"t>H_{0}")

        epi_counts = reduce(vcat,map(tr->tr.epistasis[tr.H0:end],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_uH0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_uH0)

        ###############################\

        bins = 50

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

        colgap!(plot_geno,Relative(0.01))
        rowgap!(plot_geno,Relative(0.05))
    
        colgap!(plot_wait, Relative(0.01))
        rowgap!(plot_wait, Relative(0.05))

        colgap!(plot_mut_hist, Relative(0.01))
        rowgap!(plot_mut_hist, Relative(0.05))

        # colgap!(plot_epi_types, Relative(0.01))
        # rowgap!(plot_epi_types, Relative(0.05))
    end

    if wait_time_summary == :mean
        labels_wait =  [L"\mathbb{E}[\text{Average total generations}]"]
    else
        labels_wait =  [L"\text{Average total generations}"]
    end

    labels_mut =  [L"\text{new/+}",L"\text{existing/+}",L"\text{new/*}",L"\text{existing/*}"]

    labels_epi  = [L"\text{R.S.E}",L"\text{S.E}",L"\text{N.E}",L"\text{S.M}"]

    labels = vcat(labels_wait,labels_epi)

    symbol_wait = [[wt_s_list[1], wt_l_list[1]]]

    symbol_mut = [PolyElement(color=c) for c in palette(:viridis, 4)[1:4]]

    symbol_epi = [PolyElement(color=c) for c in evo_config.pie_colors]

    symbol_all = vcat(symbol_wait,symbol_epi)

    # Legend(fig[top_n+1, :], symbol_all, labels, framevisible=false,nbanks = 1,orientation = :horizontal,patchsize = (10, 10), rowgap = 10,colgap = 10)

    # Legend(fig[top_n, 4:13, Bottom()], symbol_all, labels, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = 10,colgap = 10)

    legend_row_gap = 2

    # Legend(fig[top_n, 2:4, Bottom()], symbol_wait, labels_wait, framevisible=false,nbanks =1,orientation = :horizontal,patchsize = (10, 10),rowgap = legend_row_gap,colgap = 2,padding=(10.0f0, 10.0f0, 0f0, evo_config.fontsize+1.5*legend_row_gap))

    Legend(fig[top_n, 2:4, Bottom()], symbol_all, labels, framevisible=false,nbanks =3,orientation = :horizontal,patchsize = (10, 10),rowgap = legend_row_gap,colgap = 2,padding=(10.0f0, 10.0f0, 0f0, evo_config.fontsize+1.5*legend_row_gap))

    Legend(fig[top_n, 5:9, Bottom()], symbol_mut, labels_mut, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = legend_row_gap,colgap = 2)

    # Legend(fig[top_n,  10:13, Bottom()], symbol_epi, labels_epi, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = legend_row_gap,colgap = 2)

    linkyaxes!(ax_wait_list...)
    linkyaxes!(ax_wait_2_list...)

    rowgap!(fig.layout, Relative(0.01))
    colgap!(fig.layout, Relative(0.01))

end

function create_epi_summary_portrait!(fig,trajectories,top_n,mutation_operator::Union{MutationOperatorDual,MutationOperatorUniform},sorted_uep, vertex_top_map,wait_time_summary,evo_config)

    all_wait_times = reduce(hcat,[average_wait_time(tr) for tr in trajectories]);

    ax_wait_list = []

    ax_wait_list = []
    ax_wait_2_list = []

    wt_l_list = []
    wt_s_list = []

    min_t_u = -mutation_operator.max_w
    max_t_u = mutation_operator.max_w

    grid_values = Tuple.(findall(ones(Int(floor(sqrt(top_n))),Int(floor(sqrt(top_n)))) .> 0))

    for n in 1:top_n

        sub_plot = fig[grid_values[n]...] = GridLayout()

        if n==1
            ax_geno = Axis(sub_plot[1,1],backgroundcolor = (evo_config.color_scheme[n],evo_config.color_fade),title =L"M^{(i)}_{N_i}",aspect = DataAspect())
        else
            ax_geno = Axis(sub_plot[1,1],backgroundcolor = (evo_config.color_scheme[n],evo_config.color_fade),aspect = DataAspect())
        end

        ax_geno_s = Axis(sub_plot[2,1],backgroundcolor = :transparent,
        leftspinevisible = false,
        rightspinevisible = false,
        bottomspinevisible = false,
        topspinevisible = false,
        xticklabelsvisible = false, 
        yticklabelsvisible = false,
        xgridcolor = :transparent,
        ygridcolor = :transparent,
        xminorticksvisible = false,
        yminorticksvisible = false,
        xticksvisible = false,
        yticksvisible = false,
        xautolimitmargin = (0.0,0.0),
        yautolimitmargin = (0.0,0.0))

        hidedecorations!(ax_geno_s)

        draw_grn!(ax_geno,vertex_top_map[sorted_uep[n]],evo_config.draw_config,evo_config.node_colors,evo_config.fontsize,false,false)

        ax_wait = Axis(sub_plot[2,2:4],yticklabelsize = 0.8*evo_config.fontsize,yaxisposition = :right, xticklabelsvisible = false, yticksize= 0.25*evo_config.fontsize)

        ax_wait_2 = Axis(sub_plot[2,2:4], yticklabelcolor = :red,yscale = log10,yticklabelsize = 0.8*evo_config.fontsize,yticksize= 0.25*evo_config.fontsize)

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

        ax_epi_lH0 = Axis(sub_plot[1,2],title = L"t<H_{0}")

        epi_counts = reduce(vcat,map(tr->tr.epistasis[1:tr.H0-2],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_lH0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_lH0)

        ax_epi_H0 = Axis(sub_plot[1,3],title = L"t=H_{0}")
 
        epi_counts = map(tr->tr.epistasis[tr.H0-1],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_H0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_H0)

        ax_epi_uH0 = Axis(sub_plot[1,4],title = L"t>H_{0}")

        epi_counts = reduce(vcat,map(tr->tr.epistasis[tr.H0:end],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_uH0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_uH0)

        ###############################\

        colgap!(sub_plot, Relative(0.01))
        rowgap!(sub_plot, Relative(0.05))

    end

    if wait_time_summary == :mean
        labels_wait =  [L"\mathbb{E}[\text{Total generations}]"]
    else
        labels_wait =  [L"\text{Total generations}"]
    end

    labels_mut =  [L"\text{new:+}",L"\text{existing:+}",L"\text{new:} \times",L"\text{existing:} \times"]

    labels_epi  = [L"\text{TD}",L"\text{SD}",L"\text{TI}",L"\text{SIC}"]

    labels = reduce(vcat,[labels_wait,labels_epi,labels_mut])

    symbol_wait = [[wt_s_list[1], wt_l_list[1]]]

    symbol_mut = [PolyElement(color=c) for c in palette(:viridis, 4)[1:4]]

    symbol_epi = [PolyElement(color=c) for c in evo_config.pie_colors]

    symbol_all = reduce(vcat,[symbol_wait,symbol_epi,symbol_mut])

    # Legend(fig[top_n+1, :], symbol_all, labels, framevisible=false,nbanks = 1,orientation = :horizontal,patchsize = (10, 10), rowgap = 10,colgap = 10)

    # Legend(fig[top_n, 4:13, Bottom()], symbol_all, labels, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = 10,colgap = 10)

    legend_row_gap = 2

    # Legend(fig[top_n, 2:4, Bottom()], symbol_wait, labels_wait, framevisible=false,nbanks =1,orientation = :horizontal,patchsize = (10, 10),rowgap = legend_row_gap,colgap = 2,padding=(10.0f0, 10.0f0, 0f0, evo_config.fontsize+1.5*legend_row_gap))

    Legend(fig[2, 1:2, Bottom()], symbol_all, labels, framevisible=false,nbanks =2,orientation = :horizontal,patchsize = (10, 10),rowgap = legend_row_gap,colgap = 2,padding=(10.0f0, 10.0f0, 0f0, evo_config.fontsize+1.5*legend_row_gap))

    # Legend(fig[2, 2, Bottom()], symbol_mut, labels_mut, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = legend_row_gap,colgap = 2)

    # Legend(fig[top_n,  10:13, Bottom()], symbol_epi, labels_epi, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = legend_row_gap,colgap = 2)

    linkyaxes!(ax_wait_list...)
    linkyaxes!(ax_wait_2_list...)

    rowgap!(fig.layout, Relative(0.01))
    colgap!(fig.layout, Relative(0.01))

end

function create_epi_summary_portrait_v1!(fig,trajectories,top_n,mutation_operator::Union{MutationOperatorDual,MutationOperatorUniform},sorted_uep, vertex_top_map,wait_time_summary,evo_config)

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

        if n==1
            ax_epi_lH0 = Axis(epi_subplot[n,1],title = L"t<S_{0}")
        else
            ax_epi_lH0 = Axis(epi_subplot[n,1])
        end

        epi_counts = reduce(vcat,map(tr->tr.epistasis[1:tr.H0-2],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_lH0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_lH0)

        if n==1
            ax_epi_H0 = Axis(epi_subplot[n,2],title = L"t=S_{0}")
        else
            ax_epi_H0 = Axis(epi_subplot[n,2])
        end
 
        epi_counts = map(tr->tr.epistasis[tr.H0-1],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_H0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_H0)

        if n==1
            ax_epi_uH0 = Axis(epi_subplot[n,3],title = L"t>S_{0}")
        else
            ax_epi_uH0 = Axis(epi_subplot[n,3])
        end

        epi_counts = reduce(vcat,map(tr->tr.epistasis[tr.H0:end],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_uH0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_uH0)

        ###############################\

    end

    # colgap!(geno_subplot, Relative(0.8))
    # rowgap!(geno_subplot, Relative(0.05))

    colgap!(wait_subplot, Relative(0.01))
    rowgap!(wait_subplot, Relative(0.05))

    colgap!(epi_subplot, Relative(0.01))
    rowgap!(epi_subplot, Relative(0.05))

    ylabelwl = Label(wait_subplot[1:top_n,1,Left()], L"\text{Average wait time}", rotation = pi/2, padding = (5.,30.,0.,0.), color = :red)
    ylabelwr = Label(wait_subplot[1:top_n,end,Right()], L"\text{Mutant composition}", rotation = pi/2, padding = (30.,2.,0.,0.))

    # title_wt = Label(wait_subplot[1,1:end,TopLeft()], L"\text{Title}", padding = (0.,0.,10.,5.))

    if wait_time_summary == :mean
        labels_wait =  [L"\mathbb{E}[\text{Total generations}]"]
    else
        labels_wait =  [L"\text{Total generations}"]
    end

    labels_mut =  [L"\text{new:+}",L"\text{existing:+}",L"\text{new:} \times",L"\text{existing:} \times"]

    labels_epi  = [L"\text{TD}",L"\text{SD}",L"\text{TI}",L"\text{SIC}"]

    labels = reduce(vcat,[labels_wait,labels_epi,labels_mut])

    symbol_wait = [[wt_s_list[1], wt_l_list[1]]]

    symbol_mut = [PolyElement(color=c) for c in palette(:viridis, 4)[1:4]]

    symbol_epi = [PolyElement(color=c) for c in evo_config.pie_colors]

    symbol_wait_all = reduce(vcat,[symbol_wait,symbol_mut])
    label_wait_all = reduce(vcat,[labels_wait,labels_mut])

    legend_row_gap = 2

    Legend(epi_subplot[top_n, :, Bottom()], symbol_epi, labels_epi, framevisible=false,nbanks = 1,orientation = :horizontal,patchsize = (10, 10), colgap = 4, padding=(0.,0.,0f0, evo_config.fontsize+1.5*legend_row_gap))

    # Legend(wait_subplot[top_n, :, Bottom()], symbol_wait_all, label_wait_all, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = 10,colgap = 10,padding=(10.0f0, 10.0f0, 0f0, evo_config.fontsize+1.5*legend_row_gap))

    Legend(wait_subplot[top_n, :, Bottom()], symbol_mut, labels_mut, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), colgap = 4, rowgap = 4, padding=(0.,0.,0f0, evo_config.fontsize+1.5*legend_row_gap))

    # # Legend(fig[top_n+1, :], symbol_all, labels, framevisible=false,nbanks = 1,orientation = :horizontal,patchsize = (10, 10), rowgap = 10,colgap = 10)

    # # Legend(fig[top_n, 4:13, Bottom()], symbol_all, labels, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = 10,colgap = 10)

    # legend_row_gap = 2

    # # Legend(fig[top_n, 2:4, Bottom()], symbol_wait, labels_wait, framevisible=false,nbanks =1,orientation = :horizontal,patchsize = (10, 10),rowgap = legend_row_gap,colgap = 2,padding=(10.0f0, 10.0f0, 0f0, evo_config.fontsize+1.5*legend_row_gap))

    # Legend(fig[2, 1:2, Bottom()], symbol_all, labels, framevisible=false,nbanks =2,orientation = :horizontal,patchsize = (10, 10),rowgap = legend_row_gap,colgap = 2,padding=(10.0f0, 10.0f0, 0f0, evo_config.fontsize+1.5*legend_row_gap))

    # # Legend(fig[2, 2, Bottom()], symbol_mut, labels_mut, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = legend_row_gap,colgap = 2)

    # # Legend(fig[top_n,  10:13, Bottom()], symbol_epi, labels_epi, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = legend_row_gap,colgap = 2)

    linkyaxes!(ax_wait_list...)
    linkyaxes!(ax_wait_2_list...)

    rowgap!(fig.layout, Relative(0.1))
    colgap!(fig.layout, Relative(0.03))

end


function create_epi_summary_portrait_v2!(fig,trajectories,top_n,mutation_operator::Union{MutationOperatorDual,MutationOperatorUniform},sorted_uep, vertex_top_map,wait_time_summary,evo_config)

    all_wait_times = reduce(hcat,[average_wait_time_ext(tr) for tr in trajectories]);

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

    for n in 1:top_n

        sub_plot = fig[n,:] = GridLayout()

        if n==1
            ax_geno = Axis(geno_subplot[n,1],backgroundcolor = (evo_config.color_scheme[n],evo_config.color_fade),title =L"M^{(i)}_{N_i}",aspect = DataAspect())
        else
            ax_geno = Axis(geno_subplot[n,1],backgroundcolor = (evo_config.color_scheme[n],evo_config.color_fade),aspect = DataAspect())
        end

        draw_grn!(ax_geno,vertex_top_map[sorted_uep[n]],evo_config.draw_config,evo_config.node_colors,evo_config.fontsize,false,false)

        if n ==1
            ax_wait = Axis(wait_subplot[n,:],yticklabelsize = 0.8*evo_config.fontsize,yaxisposition = :right, xticklabelsvisible = false, yticksize= 0.25*evo_config.fontsize, title = L"t=1 \text{ } t<S_{0} \text{ }  t=S_{0}  \text{ } t=S_{0}",titlesize = 10.)
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

        mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,1,1), [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

        mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

        push!(mut_type_prop_all,mut_type_prop_av)
        push!(mut_type_labels, [1,2,3,4])
        push!(mut_type_time_labels,[1,1,1,1])

        mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,2,tr.H0-2), [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n]) & (tr.H0-2 > 2),trajectories))

        mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

        push!(mut_type_prop_all,mut_type_prop_av)
        push!(mut_type_labels, [1,2,3,4])
        push!(mut_type_time_labels,[2,2,2,2])

        # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0-1,tr.H0-1), [:existing,:new,:del]),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n] ,trajectories_p))
        mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0-1,tr.H0-1), [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n] ,trajectories))

        mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

        push!(mut_type_prop_all,mut_type_prop_av)
        push!(mut_type_labels, [1,2,3,4])
        push!(mut_type_time_labels,[3,3,3,3])

        # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0,length(tr.topologies)-1), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n])  & (tr.H0 < length(tr.topologies)),trajectories_p))

        mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0,length(tr.topologies)-1), [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n])  & (tr.H0 < length(tr.topologies)),trajectories))

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

        push!(ax_wait_list,ax_wait)

        ############################ noise_distribut

        sample_id = findall(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)

        if wait_time_summary == :mean

            mean_wait = mean(all_wait_times[:,sample_id],dims = 2)[:,1]

            std_error_wait = std(all_wait_times[:,sample_id],dims = 2)[:,1] ./ sqrt(length(sample_id))

            mean_wait_type_labels = [1,2,3,4]

            wt_l = CairoMakie.lines!(ax_wait_2,mean_wait_type_labels,mean_wait,color = :red,linewidth = evo_config.wait_linewidth)
            wt_s = CairoMakie.scatter!(ax_wait_2,mean_wait_type_labels,mean_wait,color = :red,markersize = evo_config.wait_markersize)

            CairoMakie.errorbars!(ax_wait_2,1:length(mean_wait),mean_wait,5 * std_error_wait,color = :red,whiskerwidth = evo_config.wait_markersize/2)

        else

            median_wait_time = mapslices(row->quantile(row, [0.5]),all_wait_times[:,sample_id],dims =2)[:,1]
            lq_wait_time = mapslices(row->quantile(row, [0.25]),all_wait_times[:,sample_id],dims =2)[:,1]
            uq_wait_time = mapslices(row->quantile(row, [0.75]),all_wait_times[:,sample_id],dims =2)[:,1]

            median_wait_type_labels = [1,2,3,4]

            wt_l = CairoMakie.lines!(ax_wait_2,median_wait_type_labels,median_wait_time,color = :red,linewidth = evo_config.wait_linewidth)
            wt_s = CairoMakie.scatter!(ax_wait_2,median_wait_type_labels,median_wait_time,color = :red,markersize = evo_config.wait_markersize)

            CairoMakie.rangebars!(ax_wait_2,1:length(median_wait_time),lq_wait_time,uq_wait_time,color = :red,whiskerwidth = evo_config.wait_markersize/2)
        end

        push!(ax_wait_2_list,ax_wait_2)

        push!(wt_l_list,wt_l)
        push!(wt_s_list,wt_s)

        #############################

        if n==1
            ax_epi_fH0 = Axis(epi_subplot[n,1],title = L"t=1",titlesize = 10.)
        else
            ax_epi_fH0 = Axis(epi_subplot[n,1])
        end

        epi_counts = reduce(vcat,map(tr->tr.epistasis[1],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_fH0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_fH0)

        if n==1
            ax_epi_lH0 = Axis(epi_subplot[n,2],title = L"1<t<S_{0}",titlesize = 10.)
        else
            ax_epi_lH0 = Axis(epi_subplot[n,2])
        end

        epi_counts = reduce(vcat,map(tr->tr.epistasis[2:tr.H0-2],filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n]) & (tr.H0-2 > 2),trajectories)))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_lH0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_lH0)

        if n==1
            ax_epi_H0 = Axis(epi_subplot[n,3],title = L"t=S_{0}",titlesize = 10.)
        else
            ax_epi_H0 = Axis(epi_subplot[n,3])
        end
 
        epi_counts = map(tr->tr.epistasis[tr.H0-1],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_H0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_H0)

        if n==1
            ax_epi_uH0 = Axis(epi_subplot[n,4],title = L"t>S_{0}",titlesize = 10.)
        else
            ax_epi_uH0 = Axis(epi_subplot[n,4])
        end

        epi_counts = reduce(vcat,map(tr->tr.epistasis[tr.H0:end],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_uH0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_uH0)

        ###############################\

    end

    # colgap!(geno_subplot, Relative(0.8))
    # rowgap!(geno_subplot, Relative(0.05))

    colgap!(wait_subplot, Relative(0.01))
    rowgap!(wait_subplot, Relative(0.05))

    colgap!(epi_subplot, Relative(0.01))
    rowgap!(epi_subplot, Relative(0.05))

    ylabelwl = Label(wait_subplot[1:top_n,1,Left()], L"\text{Average wait time}", rotation = pi/2, padding = (5.,30.,0.,0.), color = :red)
    ylabelwr = Label(wait_subplot[1:top_n,end,Right()], L"\text{Mutant composition}", rotation = pi/2, padding = (30.,2.,0.,0.))

    # title_wt = Label(wait_subplot[1,1:end,TopLeft()], L"\text{Title}", padding = (0.,0.,10.,5.))

    if wait_time_summary == :mean
        labels_wait =  [L"\mathbb{E}[\text{Total generations}]"]
    else
        labels_wait =  [L"\text{Total generations}"]
    end

    labels_mut =  [L"\text{new:+}",L"\text{existing:+}",L"\text{new:} \times",L"\text{existing:} \times"]

    labels_epi  = [L"\text{TD}",L"\text{SD}",L"\text{TI}",L"\text{SIC}"]

    labels = reduce(vcat,[labels_wait,labels_epi,labels_mut])

    symbol_wait = [[wt_s_list[1], wt_l_list[1]]]

    symbol_mut = [PolyElement(color=c) for c in palette(:viridis, 4)[1:4]]

    symbol_epi = [PolyElement(color=c) for c in evo_config.pie_colors]

    symbol_wait_all = reduce(vcat,[symbol_wait,symbol_mut])
    label_wait_all = reduce(vcat,[labels_wait,labels_mut])

    legend_row_gap = 2

    Legend(epi_subplot[top_n, :, Bottom()], symbol_epi, labels_epi, framevisible=false,nbanks = 1,orientation = :horizontal,patchsize = (10, 10), colgap = 4, padding=(0.,0.,0f0, evo_config.fontsize+1.5*legend_row_gap))

    # Legend(wait_subplot[top_n, :, Bottom()], symbol_wait_all, label_wait_all, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = 10,colgap = 10,padding=(10.0f0, 10.0f0, 0f0, evo_config.fontsize+1.5*legend_row_gap))

    Legend(wait_subplot[top_n, :, Bottom()], symbol_mut, labels_mut, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), colgap = 4, rowgap = 4, padding=(0.,0.,0f0, evo_config.fontsize+1.5*legend_row_gap))

    # # Legend(fig[top_n+1, :], symbol_all, labels, framevisible=false,nbanks = 1,orientation = :horizontal,patchsize = (10, 10), rowgap = 10,colgap = 10)

    # # Legend(fig[top_n, 4:13, Bottom()], symbol_all, labels, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = 10,colgap = 10)

    # legend_row_gap = 2

    # # Legend(fig[top_n, 2:4, Bottom()], symbol_wait, labels_wait, framevisible=false,nbanks =1,orientation = :horizontal,patchsize = (10, 10),rowgap = legend_row_gap,colgap = 2,padding=(10.0f0, 10.0f0, 0f0, evo_config.fontsize+1.5*legend_row_gap))

    # Legend(fig[2, 1:2, Bottom()], symbol_all, labels, framevisible=false,nbanks =2,orientation = :horizontal,patchsize = (10, 10),rowgap = legend_row_gap,colgap = 2,padding=(10.0f0, 10.0f0, 0f0, evo_config.fontsize+1.5*legend_row_gap))

    # # Legend(fig[2, 2, Bottom()], symbol_mut, labels_mut, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = legend_row_gap,colgap = 2)

    # # Legend(fig[top_n,  10:13, Bottom()], symbol_epi, labels_epi, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = legend_row_gap,colgap = 2)

    linkyaxes!(ax_wait_list...)
    linkyaxes!(ax_wait_2_list...)

    rowgap!(fig.layout, Relative(0.1))
    colgap!(fig.layout, Relative(0.03))

end

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

function create_epi_summary_portrait_bar!(fig,trajectories,top_n,mutation_operator::Union{MutationOperatorDual,MutationOperatorUniform},mut_prob,sorted_uep, vertex_top_map,wait_time_summary,evo_config)

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

function create_epi_single_portrait_bar!(fig,trajectories,mutation_operator::Union{MutationOperatorDual,MutationOperatorUniform},mut_prob,sorted_uep, vertex_top_map,wait_time_summary,evo_config)

    all_wait_times = reduce(hcat,[average_wait_time(tr) for tr in trajectories]);

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

    null_mut_n = [prob_k_mutations(k,mut_prob,10) for k in 1:10];

    Label

    ax_epi_list = []

    ax_wait = Axis(wait_subplot[1,1],yticklabelsize = 0.8*evo_config.fontsize,yaxisposition = :right, xticklabelsvisible = false, yticksize= 0.25*evo_config.fontsize, title = L"t<S_{0} \text{  }  t=S_{0}  \text{  } t=S_{0}")

    ax_wait_2 = Axis(wait_subplot[1,1], yticklabelcolor = :red,yscale = log10,yticklabelsize = 0.8*evo_config.fontsize,yticksize= 0.25*evo_config.fontsize)

    hidespines!(ax_wait_2 )
    hideydecorations!(ax_wait_2,label = false,ticklabels = false,ticks = false,minorticks = false)
    hidexdecorations!(ax_wait_2)

    # #############################

    mut_type_prop_all = []
    mut_type_time_labels = []
    mut_type_labels = []

    # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,1,tr.H0-2), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n]) & (tr.H0-2 > 0),trajectories_p))

    mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,1,tr.H0-2), [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]),filter(tr->(tr.H0-2 > 0),trajectories))

    mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

    push!(mut_type_prop_all,mut_type_prop_av)
    push!(mut_type_labels, [1,2,3,4])
    push!(mut_type_time_labels,[1,1,1,1])

    # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0-1,tr.H0-1), [:existing,:new,:del]),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n] ,trajectories_p))
    mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0-1,tr.H0-1), [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]),trajectories)

    mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

    push!(mut_type_prop_all,mut_type_prop_av)
    push!(mut_type_labels, [1,2,3,4])
    push!(mut_type_time_labels,[2,2,2,2])

    # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0,length(tr.topologies)-1), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n])  & (tr.H0 < length(tr.topologies)),trajectories_p))

    mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0,length(tr.topologies)-1), [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]),filter(tr->(tr.H0 < length(tr.topologies)),trajectories))

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

    ############################ noise_distribut

    if wait_time_summary == :mean

        mean_wait = mean(all_wait_times,dims = 2)[:,1]

        std_error_wait = std(all_wait_times,dims = 2)[:,1] ./ sqrt(length(sample_id))

        mean_wait_type_labels = [1,2,3]

        wt_l = CairoMakie.lines!(ax_wait_2,mean_wait_type_labels,mean_wait,color = :red,linewidth = evo_config.wait_linewidth)
        wt_s = CairoMakie.scatter!(ax_wait_2,mean_wait_type_labels,mean_wait,color = :red,markersize = evo_config.wait_markersize)

        CairoMakie.errorbars!(ax_wait_2,1:length(mean_wait),mean_wait,5 * std_error_wait,color = :red,whiskerwidth = evo_config.wait_markersize/2)

    else

        median_wait_time = mapslices(row->quantile(row, [0.5]),all_wait_times,dims =2)[:,1]
        lq_wait_time = mapslices(row->quantile(row, [0.25]),all_wait_times,dims =2)[:,1]
        uq_wait_time = mapslices(row->quantile(row, [0.75]),all_wait_times,dims =2)[:,1]

        median_wait_type_labels = [1,2,3]

        wt_l = CairoMakie.lines!(ax_wait_2,median_wait_type_labels,median_wait_time,color = :red,linewidth = evo_config.wait_linewidth)
        wt_s = CairoMakie.scatter!(ax_wait_2,median_wait_type_labels,median_wait_time,color = :red,markersize = evo_config.wait_markersize)

        CairoMakie.rangebars!(ax_wait_2,1:length(median_wait_time),lq_wait_time,uq_wait_time,color = :red,whiskerwidth = evo_config.wait_markersize/2)
    end

    push!(wt_l_list,wt_l)
    push!(wt_s_list,wt_s)

    #############################

    ax_epi = Axis(epi_subplot[1,1],yticklabelsize = 10.,xticklabelsize = 10,xlabel = L"\text{weight changes per mutant}",ygridvisible = false,xgridvisible = false)

    epi_counts_lS0 = reduce(vcat,map(tr->tr.epistasis[1:tr.H0-2],filter(tr->(tr.H0-2 > 0),trajectories)));
    mut_n_counts_lS0 = reduce(vcat,map(tr->get_mut_n(tr,1,tr.H0-2),filter(tr->(tr.H0-2 > 0),trajectories)));

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

    epi_counts_S0 = reduce(vcat,map(tr->tr.epistasis[tr.H0-1],trajectories));
    mut_n_counts_S0 = reduce(vcat,map(tr->get_mut_n(tr,tr.H0-1,tr.H0-1),trajectories));

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

    epi_counts_hS0 = reduce(vcat,map(tr->tr.epistasis[tr.H0:end],trajectories));
    mut_n_counts_hS0 = reduce(vcat,map(tr->get_mut_n(tr,tr.H0,length(tr.geno_traj)-1),trajectories));

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

    scatter_epi_mutn_text = reduce(vcat,[[L"n < S_0" for i in 1:10],[L"n = S_0" for i in 1:10],[L"n > S_0" for i in 1:10]]);

    CairoMakie.barplot!(ax_epi,x_epi_mutn_all,values_epi_mutn_all, color = color_epi_mutn_all, stack = stack_epi_mutn_all, dodge = dodge_epi_mutn_all,width = width)

    # rangebars!(ax_epi,xerr,first.(err_values_epi_mutn),last.(err_values_epi_mutn ); whiskerwidth = 8)

    CairoMakie.lines!(ax_epi,null_mut_n, color = :black, linestyle = :dash, linewidth = evo_config.wait_linewidth)
    CairoMakie.scatter!(ax_epi,null_mut_n, color = :black,marker = 'x',markersize = evo_config.wait_markersize)

    # CairoMakie.scatter!(ax_epi,xerr, max.(1.4 .* scatter_height_epi_mutn,0.1),color = scatter_epi_mutn_color,marker = scatter_epi_mutn_mark,markersize = evo_config.wait_markersize)

    CairoMakie.text!(ax_epi,xerr, max.(1.1 .* scatter_height_epi_mutn,0.1),text = scatter_epi_mutn_text,rotation = pi/2,fontsize = 10.,align = (:left, :center))

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

    ax_wait = Axis(wait_subplot[1,1],yticklabelsize = 0.8*evo_config.fontsize,yaxisposition = :right, xticklabelsvisible = false, yticksize= 0.25*evo_config.fontsize, title = L"n=1 \text{ } 1<n<S_{0} \text{ }  n=S_{0}  \text{ } n=S_{0}")

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

    epi_counts_1 = reduce(vcat,map(tr->tr.epistasis[1],filter(tr->(tr.H0 > 2),trajectories)));
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

    epi_counts_lS0 = reduce(vcat,map(tr->tr.epistasis[2:tr.H0-2],filter(tr->(tr.H0-2 > 0),trajectories)));
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

    epi_counts_S0 = reduce(vcat,map(tr->tr.epistasis[tr.H0-1],trajectories));
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

    epi_counts_hS0 = reduce(vcat,map(tr->tr.epistasis[tr.H0:end],trajectories));
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

function create_extensive_epi_summary_portrait_v1!(fig,trajectories,top_n,mutation_operator::Union{MutationOperatorDual,MutationOperatorUniform},sorted_uep, vertex_top_map,wait_time_summary,evo_config)

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

        mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,1,1), [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n] ,trajectories))

        mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

        push!(mut_type_prop_all,mut_type_prop_av)
        push!(mut_type_labels, [1,2,3,4])
        push!(mut_type_time_labels,[1,1,1,1])

        # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0-1,tr.H0-1), [:existing,:new,:del]),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n] ,trajectories_p))
        mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0-2,tr.H0-2), [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n]) & (tr.H0-2 > 0),trajectories))

        mut_type_prop_av = mean(reduce(hcat,mut_type_prop),dims = 2)[:,1]

        push!(mut_type_prop_all,mut_type_prop_av)
        push!(mut_type_labels, [1,2,3,4])
        push!(mut_type_time_labels,[2,2,2,2])

        # mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0,length(tr.topologies)-1), [:existing,:new,:del]),filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n])  & (tr.H0 < length(tr.topologies)),trajectories_p))

        mut_type_prop = map(tr->calculate_mut_type_proportion(get_mut_type(tr,tr.H0-1,tr.H0-1), [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n] ,trajectories))

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

        if n==1
            ax_epi_lH0 = Axis(epi_subplot[n,1],title = L"t=1")
        else
            ax_epi_lH0 = Axis(epi_subplot[n,1])
        end

        epi_counts = map(tr->tr.epistasis[1],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_lH0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_lH0)

        if n==1
            ax_epi_H0 = Axis(epi_subplot[n,2],title = L"t=S_{0}-1")
        else
            ax_epi_H0 = Axis(epi_subplot[n,2])
        end
 
        epi_counts = map(tr->tr.epistasis[tr.H0-2],filter(tr->(tr.inc_metagraph_vertices[end] == sorted_uep[n]) & (tr.H0-2 > 0),trajectories))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_H0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_H0)

        if n==1
            ax_epi_uH0 = Axis(epi_subplot[n,3],title = L"t=S_{0}")
        else
            ax_epi_uH0 = Axis(epi_subplot[n,3])
        end

        epi_counts = reduce(vcat,map(tr->tr.epistasis[tr.H0-1],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_uH0,epi_counts_prop,radius = evo_config.pie_radius,color = evo_config.pie_colors,
        inner_radius = evo_config.pie_inner_radius,
        strokecolor = :white,
        strokewidth = evo_config.pie_strokewidth)

        CairoMakie.hidedecorations!(ax_epi_uH0)

        ###############################\

    end

    # colgap!(geno_subplot, Relative(0.8))
    # rowgap!(geno_subplot, Relative(0.05))

    colgap!(wait_subplot, Relative(0.01))
    rowgap!(wait_subplot, Relative(0.05))

    colgap!(epi_subplot, Relative(0.01))
    rowgap!(epi_subplot, Relative(0.05))

    ylabelwl = Label(wait_subplot[1:top_n,1,Left()], L"\text{Average wait time}", rotation = pi/2, padding = (5.,30.,0.,0.), color = :red)
    ylabelwr = Label(wait_subplot[1:top_n,end,Right()], L"\text{Mutant composition}", rotation = pi/2, padding = (30.,2.,0.,0.))

    # title_wt = Label(wait_subplot[1,1:end,TopLeft()], L"\text{Title}", padding = (0.,0.,10.,5.))

    if wait_time_summary == :mean
        labels_wait =  [L"\mathbb{E}[\text{Total generations}]"]
    else
        labels_wait =  [L"\text{Total generations}"]
    end

    labels_mut =  [L"\text{new:+}",L"\text{existing:+}",L"\text{new:} \times",L"\text{existing:} \times"]

    labels_epi  = [L"\text{TD}",L"\text{SD}",L"\text{TI}",L"\text{SIC}"]

    labels = reduce(vcat,[labels_wait,labels_epi,labels_mut])

    symbol_wait = [[wt_s_list[1], wt_l_list[1]]]

    symbol_mut = [PolyElement(color=c) for c in palette(:viridis, 4)[1:4]]

    symbol_epi = [PolyElement(color=c) for c in evo_config.pie_colors]

    symbol_wait_all = reduce(vcat,[symbol_wait,symbol_mut])
    label_wait_all = reduce(vcat,[labels_wait,labels_mut])

    legend_row_gap = 2

    Legend(epi_subplot[top_n, :, Bottom()], symbol_epi, labels_epi, framevisible=false,nbanks = 1,orientation = :horizontal,patchsize = (10, 10), colgap = 4, padding=(0.,0.,0f0, evo_config.fontsize+1.5*legend_row_gap))

    # Legend(wait_subplot[top_n, :, Bottom()], symbol_wait_all, label_wait_all, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = 10,colgap = 10,padding=(10.0f0, 10.0f0, 0f0, evo_config.fontsize+1.5*legend_row_gap))

    Legend(wait_subplot[top_n, :, Bottom()], symbol_mut, labels_mut, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), colgap = 4, rowgap = 4, padding=(0.,0.,0f0, evo_config.fontsize+1.5*legend_row_gap))

    # # Legend(fig[top_n+1, :], symbol_all, labels, framevisible=false,nbanks = 1,orientation = :horizontal,patchsize = (10, 10), rowgap = 10,colgap = 10)

    # # Legend(fig[top_n, 4:13, Bottom()], symbol_all, labels, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = 10,colgap = 10)

    # legend_row_gap = 2

    # # Legend(fig[top_n, 2:4, Bottom()], symbol_wait, labels_wait, framevisible=false,nbanks =1,orientation = :horizontal,patchsize = (10, 10),rowgap = legend_row_gap,colgap = 2,padding=(10.0f0, 10.0f0, 0f0, evo_config.fontsize+1.5*legend_row_gap))

    # Legend(fig[2, 1:2, Bottom()], symbol_all, labels, framevisible=false,nbanks =2,orientation = :horizontal,patchsize = (10, 10),rowgap = legend_row_gap,colgap = 2,padding=(10.0f0, 10.0f0, 0f0, evo_config.fontsize+1.5*legend_row_gap))

    # # Legend(fig[2, 2, Bottom()], symbol_mut, labels_mut, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = legend_row_gap,colgap = 2)

    # # Legend(fig[top_n,  10:13, Bottom()], symbol_epi, labels_epi, framevisible=false,nbanks = 2,orientation = :horizontal,patchsize = (10, 10), rowgap = legend_row_gap,colgap = 2)

    linkyaxes!(ax_wait_list...)
    linkyaxes!(ax_wait_2_list...)

    rowgap!(fig.layout, Relative(0.1))
    colgap!(fig.layout, Relative(0.03))

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

        epi_counts = reduce(vcat,map(tr->tr.epistasis[1:tr.H0-2],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

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

        epi_counts = map(tr->tr.epistasis[tr.H0-1],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

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

        epi_counts = reduce(vcat,map(tr->tr.epistasis[tr.H0:end],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

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

# function create_weight_edit_summary!(fig,n,trajectories,mutation_op::MutationOperator,sorted_uep, vertex_top_map, draw_config, node_colors,fontsize,color_scheme)

#     weight_names_latex = reshape([L"W_{aa}",L"W_{ab}",L"W_{ac}",L"W_{ba}",L"W_{bb}",L"W_{bc}",L"W_{ca}",L"W_{cb}",L"W_{cc}",L"W_{ma}",L"W_{mb}",L"W_{mc}"],(3,4));

#     grid_values = Tuple.(findall(ones(3,4) .> 0))

#     colors = reverse(palette(:tab10)[1:4])

#     ax_wait_list = []

#     plot_geno = fig[grid_values[11]...] = GridLayout()

#     ax_geno = Axis(plot_geno[1,1],backgroundcolor = (color_scheme[n],1.),title =L"\text{Minimal Stripe Topology}",aspect = DataAspect())

#     draw_grn!(ax_geno,vertex_top_map[sorted_uep[n]],draw_config,node_colors,fontsize,weight_names_latex,true,false)

#     norm_type = :pdf

#     ###############################

#     plot_mut_hist = fig[grid_values[12]...] = GridLayout()

#     bins = 50

#     min_t_list = [[],[]]
#     max_t_list = [[],[]]

#     weight_grid_layouts = []

#     for (nt,type) in enumerate([:new,:existing])

#         if type == :existing
#             mut_noise_dist = mutation_op.noise_distribution;
#         else
#             mut_noise_dist = Uniform(-mutation_op.max_w,mutation_op.max_w);
#         end

#         mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

#         if type == :existing
#             if n==1
#                 ax1 = Axis(plot_mut_hist[1,1],title = L"t<H_{0}")
#             else
#                 ax1 = Axis(plot_mut_hist[1,1])
#             end
#         else
#             ax1 = Axis(plot_mut_hist[2,1])
#         end

#         if type == :existing
#             CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[1])
#         else
#             CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[2])
#         end

#         min_t = minimum(mut_size)
#         max_t = maximum(mut_size)

#         push!(min_t_list[nt],min_t)
#         push!(max_t_list[nt],max_t)

#         norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];

#         CairoMakie.lines!(ax1,LinRange(min_t,max_t,100),norm_pdf,color = :red)
#         CairoMakie.vlines!(ax1,0,color = :red,linestyle = "--")

#         mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

#         if type == :existing
#             if n==1
#                 ax2 = Axis(plot_mut_hist[1,2],title = L"t=H_{0}")
#             else
#                 ax2 = Axis(plot_mut_hist[1,2])
#             end
#         else
#             ax2 = Axis(plot_mut_hist[2,2])
#         end

#         min_t = minimum(mut_size)
#         max_t = maximum(mut_size)

#         push!(min_t_list[nt],min_t)
#         push!(max_t_list[nt],max_t)

#         norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];

#         if type == :existing
#             CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[1])
#         else
#             CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[2])
#         end

#         CairoMakie.lines!(ax2,LinRange(min_t,max_t,100),norm_pdf,color = :red)
#         CairoMakie.vlines!(ax2,0,color = :red,linestyle = "--")

#         mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,tr.H0+1,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

#         if type == :existing
#             if n==1
#                 ax3 = Axis(plot_mut_hist[1,3],title = L"t>H_{0}")
#             else
#                 ax3 = Axis(plot_mut_hist[1,3])
#             end
#         else
#             ax3 = Axis(plot_mut_hist[2,3])
#         end

#         min_t = minimum(mut_size)
#         max_t = maximum(mut_size)

#         push!(min_t_list[nt],min_t)
#         push!(max_t_list[nt],max_t)

#         norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];

#         if type == :existing
#             CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[1])
#         else
#             CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[2])
#         end

#         CairoMakie.lines!(ax3,LinRange(min_t,max_t,100),norm_pdf,color = :red)
#         CairoMakie.vlines!(ax3,0,color = :red,linestyle = "--")

#         # linkxaxes!([ax1,ax2,ax3]...)
#         # linkyaxes!([ax1,ax2,ax3]...)

#         hidedecorations!(ax1)
#         hidedecorations!(ax2)
#         hidedecorations!(ax3)
#     end

#     colgap!(plot_mut_hist, 10)
#     rowgap!(plot_mut_hist, 10)

#     push!(weight_grid_layouts,plot_mut_hist)

#     ###############################\

#     bins = 15

#     for weight_id in 1:10

#         plot_mut_hist = fig[grid_values[weight_id]...] = GridLayout()

#         for (nt,type) in enumerate([:new,:existing])

#             if type == :existing
#                 mut_noise_dist = Normal(0.0,noise_cv);
#             else
#                 mut_noise_dist = Uniform(-max_w,max_w);
#             end

#             mut_size = reduce(vcat,map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

#             if type == :existing
#                 ax1 = Axis(plot_mut_hist[1,1],title = L"t<H_{0}")
#             else
#                 ax1 = Axis(plot_mut_hist[2,1])
#             end

#             min_t = min_t_list[nt][1]
#             max_t = max_t_list[nt][1]

#             if length(mut_size) > 1
#                 if type == :existing
#                     CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[1])
#                 else
#                     CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[2])
#                 end
#             end

#             norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];
#             CairoMakie.lines!(ax1,LinRange(min_t,max_t,100),norm_pdf,color = :red)
#             CairoMakie.vlines!(ax1,0,color = :red,linestyle = "--")

#             mut_size = reduce(vcat,map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

#             if type == :existing
#                 ax2 = Axis(plot_mut_hist[1,2],title = L"t=H_{0}")
#             else
#                 ax2 = Axis(plot_mut_hist[2,2])
#             end

#             min_t = min_t_list[nt][2]
#             max_t = max_t_list[nt][2]

#             if length(mut_size) > 1
#                 if type == :existing
#                     CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[1])
#                 else
#                     CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[2])
#                 end
#             end

#             norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];
#             CairoMakie.lines!(ax2,LinRange(min_t,max_t,100),norm_pdf,color = :red)
#             CairoMakie.vlines!(ax2,0,color = :red,linestyle = "--")

#             mut_size = reduce(vcat,map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,tr.H0+1,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

#             if type == :existing
#                 ax3 = Axis(plot_mut_hist[1,3],title = L"t>H_{0}")
#             else
#                 ax3 = Axis(plot_mut_hist[2,3])
#             end

#             min_t = min_t_list[nt][3]
#             max_t = max_t_list[nt][3]

#             if length(mut_size) > 1
#                 if type == :existing
#                     CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[1])
#                 else
#                     CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[2])
#                 end
#             end

#             norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];
#             CairoMakie.lines!(ax3,LinRange(min_t,max_t,100),norm_pdf,color = :red)
#             CairoMakie.vlines!(ax3,0,color = :red,linestyle = "--")

#             # linkxaxes!([ax1,ax2,ax3]...)
#             # linkyaxes!([ax1,ax2,ax3]...)

#             hidedecorations!(ax1)
#             hidedecorations!(ax2)
#             hidedecorations!(ax3)
#         end

#         colgap!(plot_mut_hist, 10)
#         rowgap!(plot_mut_hist, 10)

#         push!(weight_grid_layouts,plot_mut_hist)

#     end

#     for (label, layout) in zip(weight_names_latex, weight_grid_layouts[2:end])
#         Label(layout[1, 1, TopLeft()], label,
#             fontsize = 26,
#             font = :bold,
#             padding = (0, 5, 5, 0),
#             halign = :right)
#     end

#     Label(weight_grid_layouts[1][1, 1, TopLeft()], L"\text{All}",
#     fontsize = 26,
#     font = :bold,
#     padding = (0, 5, 5, 0),
#     halign = :right)

#     colors = palette(:viridis, 3)

# end

# function create_weight_edit_summary!(fig,n,trajectories,mutation_op::MutationOperatorDual,sorted_uep, vertex_top_map, draw_config, node_colors,fontsize,color_scheme)

#     weight_names_latex = reshape([L"W_{aa}",L"W_{ab}",L"W_{ac}",L"W_{ba}",L"W_{bb}",L"W_{bc}",L"W_{ca}",L"W_{cb}",L"W_{cc}",L"W_{ma}",L"W_{mb}",L"W_{mc}"],(3,4));

#     grid_values = Tuple.(findall(ones(3,4) .> 0))

#     colors = reverse(palette(:tab10)[1:4])

#     ax_wait_list = []

#     plot_geno = fig[grid_values[11]...] = GridLayout()

#     ax_geno = Axis(plot_geno[1,1],backgroundcolor = (color_scheme[n],1.),title =L"\text{Minimal Stripe Topology}",aspect = DataAspect())

#     draw_grn!(ax_geno,vertex_top_map[sorted_uep[n]],draw_config,node_colors,fontsize,weight_names_latex,true,false)

#     norm_type = :pdf

#     ###############################

#     plot_mut_hist = fig[grid_values[12]...] = GridLayout()

#     bins = 50

#     min_t_list = [[],[]]
#     max_t_list = [[],[]]

#     weight_grid_layouts = []

#     for (nt,type) in enumerate([:new,:existing])

#         if type == :existing
#             mut_noise_dist = mutation_op.mult_noise_distribution;
#         else
#             mut_noise_dist = mutation_op.additive_noise_distribution;
#         end

#         mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

#         if type == :existing
#             if n==1
#                 ax1 = Axis(plot_mut_hist[1,1],title = L"t<H_{0}")
#             else
#                 ax1 = Axis(plot_mut_hist[1,1])
#             end
#         else
#             ax1 = Axis(plot_mut_hist[2,1])
#         end

#         if type == :existing
#             CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[1])
#         else
#             CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[2])
#         end

#         min_t = minimum(mut_size)
#         max_t = maximum(mut_size)

#         push!(min_t_list[nt],min_t)
#         push!(max_t_list[nt],max_t)

#         norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];

#         CairoMakie.lines!(ax1,LinRange(min_t,max_t,100),norm_pdf,color = :red)
#         CairoMakie.vlines!(ax1,0,color = :red,linestyle = "--")

#         mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

#         if type == :existing
#             if n==1
#                 ax2 = Axis(plot_mut_hist[1,2],title = L"t=H_{0}")
#             else
#                 ax2 = Axis(plot_mut_hist[1,2])
#             end
#         else
#             ax2 = Axis(plot_mut_hist[2,2])
#         end

#         min_t = minimum(mut_size)
#         max_t = maximum(mut_size)

#         push!(min_t_list[nt],min_t)
#         push!(max_t_list[nt],max_t)

#         norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];

#         if type == :existing
#             CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[1])
#         else
#             CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[2])
#         end

#         CairoMakie.lines!(ax2,LinRange(min_t,max_t,100),norm_pdf,color = :red)
#         CairoMakie.vlines!(ax2,0,color = :red,linestyle = "--")

#         mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,tr.H0+1,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

#         if type == :existing
#             if n==1
#                 ax3 = Axis(plot_mut_hist[1,3],title = L"t>H_{0}")
#             else
#                 ax3 = Axis(plot_mut_hist[1,3])
#             end
#         else
#             ax3 = Axis(plot_mut_hist[2,3])
#         end

#         min_t = minimum(mut_size)
#         max_t = maximum(mut_size)

#         push!(min_t_list[nt],min_t)
#         push!(max_t_list[nt],max_t)

#         norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];

#         if type == :existing
#             CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[1])
#         else
#             CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 3)[2])
#         end

#         CairoMakie.lines!(ax3,LinRange(min_t,max_t,100),norm_pdf,color = :red)
#         CairoMakie.vlines!(ax3,0,color = :red,linestyle = "--")

#         # linkxaxes!([ax1,ax2,ax3]...)
#         # linkyaxes!([ax1,ax2,ax3]...)

#         hidedecorations!(ax1)
#         hidedecorations!(ax2)
#         hidedecorations!(ax3)
#     end

#     colgap!(plot_mut_hist, 10)
#     rowgap!(plot_mut_hist, 10)

#     push!(weight_grid_layouts,plot_mut_hist)

#     ###############################\

#     bins = 15

#     for weight_id in 1:10

#         plot_mut_hist = fig[grid_values[weight_id]...] = GridLayout()

#         for (nt,type) in enumerate([:new,:existing])

#             if type == :existing
#                 mut_noise_dist = mutation_op.mult_noise_distribution;
#             else
#                 mut_noise_dist = mutation_op.additive_noise_distribution;
#             end    

#             mut_size = reduce(vcat,map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

#             if type == :existing
#                 ax1 = Axis(plot_mut_hist[1,1],title = L"t<H_{0}")
#             else
#                 ax1 = Axis(plot_mut_hist[2,1])
#             end

#             min_t = min_t_list[nt][1]
#             max_t = max_t_list[nt][1]

#             if length(mut_size) > 1
#                 if type == :existing
#                     CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[1])
#                 else
#                     CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[2])
#                 end
#             end

#             norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];
#             CairoMakie.lines!(ax1,LinRange(min_t,max_t,100),norm_pdf,color = :red)
#             CairoMakie.vlines!(ax1,0,color = :red,linestyle = "--")

#             mut_size = reduce(vcat,map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

#             if type == :existing
#                 ax2 = Axis(plot_mut_hist[1,2],title = L"t=H_{0}")
#             else
#                 ax2 = Axis(plot_mut_hist[2,2])
#             end

#             min_t = min_t_list[nt][2]
#             max_t = max_t_list[nt][2]

#             if length(mut_size) > 1
#                 if type == :existing
#                     CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[1])
#                 else
#                     CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[2])
#                 end
#             end

#             norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];
#             CairoMakie.lines!(ax2,LinRange(min_t,max_t,100),norm_pdf,color = :red)
#             CairoMakie.vlines!(ax2,0,color = :red,linestyle = "--")

#             mut_size = reduce(vcat,map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,tr.H0+1,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

#             if type == :existing
#                 ax3 = Axis(plot_mut_hist[1,3],title = L"t>H_{0}")
#             else
#                 ax3 = Axis(plot_mut_hist[2,3])
#             end

#             min_t = min_t_list[nt][3]
#             max_t = max_t_list[nt][3]

#             if length(mut_size) > 1
#                 if type == :existing
#                     CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[1])
#                 else
#                     CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = norm_type,color = palette(:viridis, 3)[2])
#                 end
#             end

#             norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];
#             CairoMakie.lines!(ax3,LinRange(min_t,max_t,100),norm_pdf,color = :red)
#             CairoMakie.vlines!(ax3,0,color = :red,linestyle = "--")

#             # linkxaxes!([ax1,ax2,ax3]...)
#             # linkyaxes!([ax1,ax2,ax3]...)

#             hidedecorations!(ax1)
#             hidedecorations!(ax2)
#             hidedecorations!(ax3)
#         end

#         colgap!(plot_mut_hist, 10)
#         rowgap!(plot_mut_hist, 10)

#         push!(weight_grid_layouts,plot_mut_hist)

#     end

#     for (label, layout) in zip(weight_names_latex, weight_grid_layouts[2:end])
#         Label(layout[1, 1, TopLeft()], label,
#             fontsize = 26,
#             font = :bold,
#             padding = (0, 5, 5, 0),
#             halign = :right)
#     end

#     Label(weight_grid_layouts[1][1, 1, TopLeft()], L"\text{All}",
#     fontsize = 26,
#     font = :bold,
#     padding = (0, 5, 5, 0),
#     halign = :right)

#     colors = palette(:viridis, 3)

# end

# function create_weight_edit_summary!(fig,n,trajectories,mutation_operator::MutationOperatorDual,sorted_uep, vertex_top_map, evo_config,color_scheme)

#     weight_names_latex = reshape([L"W_{aa}",L"W_{ab}",L"W_{ac}",L"W_{ba}",L"W_{bb}",L"W_{bc}",L"W_{ca}",L"W_{cb}",L"W_{cc}",L"W_{ma}",L"W_{mb}",L"W_{mc}"],(3,4));

#     grid_values = Tuple.(findall(ones(3,4) .> 0))

#     colors = reverse(palette(:tab10)[1:4])

#     ax_wait_list = []

#     plot_geno = fig[grid_values[11]...] = GridLayout()

#     ax_geno = Axis(plot_geno[1,1],backgroundcolor = (color_scheme[n],1.),title =L"\text{Minimal Stripe Topology}",aspect = DataAspect())

#     draw_grn!(ax_geno,vertex_top_map[sorted_uep[n]],evo_config.draw_config,evo_config.node_colors,evo_config.fontsize,true,false)

#     norm_type = :pdf

#     ###############################

#     plot_mut_hist = fig[grid_values[12]...] = GridLayout()

#     bins = 50

#     min_t_list = [[],[]]
#     max_t_list = [[],[]]

#     weight_grid_layouts = []

#     additive_ax = []
#     mult_ax = []

#     for type in [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]

#         if type[2] == :additive

#             min_t = -quantile(mutation_operator.additive_noise_distribution,0.99)
#             max_t = quantile(mutation_operator.additive_noise_distribution,0.99)

#             if type[1] 
#                 if n==1
#                     ax1 = Axis(plot_mut_hist[1,1],title = L"t<H_{0}")
#                 else
#                     ax1 = Axis(plot_mut_hist[1,1])
#                 end

#                 mut_noise_dist = mutation_operator.additive_noise_distribution;

#                 mut_size_r = map(tr->get_mut_size_by_type(tr,type,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))
                
#                 if length(mut_size_r) != 0
#                     mut_size = reduce(vcat,mut_size_r)
#                     CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[1]) # new additive
#                 else
#                     mut_size = [0.]
#                 end

#                 hidedecorations!(ax1)
#             else

#                 ax1 = Axis(plot_mut_hist[2,1])

#                 mut_noise_dist = mutation_operator.additive_noise_distribution;

#                 mut_size_r = map(tr->get_mut_size_by_type(tr,type,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

#                 if length(mut_size_r) != 0
#                     mut_size = reduce(vcat,mut_size_r)
#                     CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[2]) # ex additive
#                 else
#                     mut_size = [0.]
#                 end

#                 hidedecorations!(ax1)
#             end

#             norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];

#             CairoMakie.lines!(ax1,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
#             CairoMakie.vlines!(ax1,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth) # new all

#             if type[1] 
#                 if n==1
#                     ax2 = Axis(plot_mut_hist[1,2],title = L"t=H_{0}")
#                 else
#                     ax2 = Axis(plot_mut_hist[1,2])
#                 end

#                 mut_size_r = map(tr->get_mut_size_by_type(tr,type,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

#                 if length(mut_size_r) != 0
#                     mut_size = reduce(vcat,mut_size_r)
#                     CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[1])  # new additive
#                 else
#                     mut_size = [0.]
#                 end

#                 hidedecorations!(ax2)
#             else
#                 ax2 = Axis(plot_mut_hist[2,2])

#                 mut_size_r = map(tr->get_mut_size_by_type(tr,type,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

#                 if length(mut_size_r) != 0
#                     mut_size = reduce(vcat,mut_size_r)
#                     CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[2]) # ex additive
#                 else
#                     mut_size = [0.]
#                 end

#                 hidedecorations!(ax2)
#             end

#             norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];

#             CairoMakie.lines!(ax2,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
#             CairoMakie.vlines!(ax2,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth) # new all

#             if type[1] 
#                 if n==1
#                     ax3 = Axis(plot_mut_hist[1,3],title = L"t>H_{0}")
#                 else
#                     ax3 = Axis(plot_mut_hist[1,3])
#                 end

#                 mut_size_r = map(tr->get_mut_size_by_type(tr,type,tr.H0,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

#                 if length(mut_size_r) != 0
#                     mut_size = reduce(vcat,mut_size_r)
#                     CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[1]) # new additive
#                 else
#                     mut_size = [0.]
#                 end

#                 hidedecorations!(ax3)
#             else
#                 ax3 = Axis(plot_mut_hist[2,3])

#                 mut_size_r = map(tr->get_mut_size_by_type(tr,type,tr.H0,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

#                 if length(mut_size_r) != 0
#                     mut_size = reduce(vcat,mut_size_r)
#                     CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[2]) # ex additive
#                 else
#                     mut_size = [0.]
#                 end

#                 hidedecorations!(ax3)
#             end

#             norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];

#             CairoMakie.lines!(ax3,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
#             CairoMakie.vlines!(ax3,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth)

#             # hidedecorations!(ax1)
#             # hidedecorations!(ax2)
#             # hidedecorations!(ax3) # new all

#             push!(additive_ax,ax1)
#             push!(additive_ax,ax2)
#             push!(additive_ax,ax3)

#         else
#             min_t = 0.
#             max_t = quantile(mutation_operator.mult_noise_distribution,0.99)

#             if type[1]
#                 ax1 = Axis(plot_mut_hist[3,1])

#                 mut_noise_dist = mutation_operator.mult_noise_distribution;

#                 mut_size_r = map(tr->get_mut_size_by_type(tr,type,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))
                
#                 mut_size = abs.(reduce(vcat,mut_size_r))

#                 if length(mut_size) != 0
#                     CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[3]) # new mult
#                 else
#                     mut_size = [0.]
#                 end

#                 hidedecorations!(ax1)
#             else
#                 ax1 = Axis(plot_mut_hist[4,1])

#                 mut_noise_dist = mutation_operator.mult_noise_distribution;

#                 mut_size_r = map(tr->get_mut_size_by_type(tr,type,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))
                
#                 mut_size = abs.(reduce(vcat,mut_size_r))

#                 if length(mut_size) != 0
#                     CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[4]) # ex mult
#                 else
#                     mut_size = [0.]
#                 end

#                 hidedecorations!(ax1)
#             end

#             norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];

#             CairoMakie.lines!(ax1,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
#             CairoMakie.vlines!(ax1,1.,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth) # ex all

#             if type[1]
#                 ax2 = Axis(plot_mut_hist[3,2])

#                 mut_size_r = map(tr->get_mut_size_by_type(tr,type,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

#                 mut_size = abs.(reduce(vcat,mut_size_r))

#                 if length(mut_size) != 0
#                     CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[3]) # new mult
#                 else
#                     mut_size = [0.]
#                 end
                
#                 hidedecorations!(ax2)
#             else
#                 ax2 = Axis(plot_mut_hist[4,2])

#                 mut_size_r = map(tr->get_mut_size_by_type(tr,type,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

#                 mut_size = abs.(reduce(vcat,mut_size_r))

#                 if length(mut_size) != 0
#                     CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[4]) # ex mult
#                 else
#                     mut_size = [0.]
#                 end

#                 hidedecorations!(ax2)
#             end

#             norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];

#             CairoMakie.lines!(ax2,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
#             CairoMakie.vlines!(ax2,1.,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth) # ex all

#             if type[1]
#                 ax3 = Axis(plot_mut_hist[3,3])

#                 mut_size_r = map(tr->get_mut_size_by_type(tr,type,tr.H0,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

#                 mut_size = abs.(reduce(vcat,mut_size_r))

#                 if length(mut_size) != 0
#                     CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[3]) # new mult
#                 else
#                     mut_size = [0.]
#                 end

#                 hidedecorations!(ax3) # ex all
#             else
#                 ax3 = Axis(plot_mut_hist[4,3])

#                 mut_size_r = map(tr->get_mut_size_by_type(tr,type,tr.H0,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

#                 mut_size = abs.(reduce(vcat,mut_size_r))

#                 if length(mut_size) != 0
#                     CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[4]) # ex mult
#                 else
#                     mut_size = [0.]
#                 end
                
#                 hidedecorations!(ax3) # ex all
#             end

#             norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];

#             CairoMakie.lines!(ax3,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
#             CairoMakie.vlines!(ax3,1.,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth)

#             push!(mult_ax,ax1)
#             push!(mult_ax,ax2)
#             push!(mult_ax,ax3)
#         end

#     end

#     linkxaxes!(additive_ax...)
#     linkxaxes!(mult_ax...)

#     colgap!(plot_mut_hist, 10)
#     rowgap!(plot_mut_hist, 10)

#     push!(weight_grid_layouts,plot_mut_hist)

#     ###############################\ get_mut_size_by_type_and_weight(tr,type,weight_id,1,tr.H0-2)

#     bins = 15

#     for weight_id in 1:10

#         plot_mut_hist = fig[grid_values[weight_id]...] = GridLayout()

#         additive_ax = []
#         mult_ax = []
    
#         for type in [(true,:additive),(false,:additive),(true,:multiplicative),(false,:multiplicative)]
    
#             if type[2] == :additive

#                 min_t = -quantile(mutation_operator.additive_noise_distribution,0.99)
#                 max_t = quantile(mutation_operator.additive_noise_distribution,0.99)
    
#                 if type[1] 
#                     if n==1
#                         ax1 = Axis(plot_mut_hist[1,1],title = L"t<H_{0}")
#                     else
#                         ax1 = Axis(plot_mut_hist[1,1])
#                     end
    
#                     mut_noise_dist = mutation_operator.additive_noise_distribution;
    
#                     mut_size_r = map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))
                    
#                     mut_size = reduce(vcat,mut_size_r)
                    
#                     if length(mut_size) != 0
#                         CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[1]) # new additive
#                     else
#                         mut_size = [0.]
#                     end
    
#                     hidedecorations!(ax1)
#                 else
    
#                     ax1 = Axis(plot_mut_hist[2,1])
    
#                     mut_noise_dist = mutation_operator.additive_noise_distribution;
    
#                     mut_size_r = map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))
    
#                     mut_size = reduce(vcat,mut_size_r)
                    
#                     if length(mut_size) != 0
#                         CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[2]) # ex additive
#                     else
#                         mut_size = [0.]
#                     end
    
#                     hidedecorations!(ax1)
#                 end
    
#                 norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];
    
#                 CairoMakie.lines!(ax1,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
#                 CairoMakie.vlines!(ax1,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth) # new all
    
#                 if type[1] 
#                     if n==1
#                         ax2 = Axis(plot_mut_hist[1,2],title = L"t=H_{0}")
#                     else
#                         ax2 = Axis(plot_mut_hist[1,2])
#                     end
    
#                     mut_size_r = map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))
    
#                     mut_size = reduce(vcat,mut_size_r)
                    
#                     if length(mut_size) != 0
#                         CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[1])  # new additive
#                     else
#                         mut_size = [0.]
#                     end
    
#                     hidedecorations!(ax2)
#                 else
#                     ax2 = Axis(plot_mut_hist[2,2])
    
#                     mut_size_r = map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))
    
#                     mut_size = reduce(vcat,mut_size_r)
                    
#                     if length(mut_size) != 0
#                         CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[2]) # ex additive
#                     else
#                         mut_size = [0.]
#                     end
    
#                     hidedecorations!(ax2)
#                 end
    
#                 norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];
    
#                 CairoMakie.lines!(ax2,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
#                 CairoMakie.vlines!(ax2,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth) # new all
    
#                 if type[1] 
#                     if n==1
#                         ax3 = Axis(plot_mut_hist[1,3],title = L"t>H_{0}")
#                     else
#                         ax3 = Axis(plot_mut_hist[1,3])
#                     end
    
#                     mut_size_r = map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,tr.H0,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))
    
#                     mut_size = reduce(vcat,mut_size_r)
                    
#                     if length(mut_size) != 0

#                         CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[1]) # new additive
#                     else
#                         mut_size = [0.]
#                     end
    
#                     hidedecorations!(ax3)
#                 else
#                     ax3 = Axis(plot_mut_hist[2,3])
    
#                     mut_size_r = map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,tr.H0,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))
    
#                     mut_size = reduce(vcat,mut_size_r)

#                     if length(mut_size) != 0
#                         CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[2]) # ex additive
#                     else
#                         mut_size = [0.]
#                     end
    
#                     hidedecorations!(ax3)
#                 end
    
#                 norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];
    
#                 CairoMakie.lines!(ax3,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
#                 CairoMakie.vlines!(ax3,0,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth)
    
#                 # hidedecorations!(ax1)
#                 # hidedecorations!(ax2)
#                 # hidedecorations!(ax3) # new all
    
#                 push!(additive_ax,ax1)
#                 push!(additive_ax,ax2)
#                 push!(additive_ax,ax3)
    
#             else
#                 min_t = 0.
#                 max_t = quantile(mutation_operator.mult_noise_distribution,0.99)
    
#                 if type[1]
#                     ax1 = Axis(plot_mut_hist[3,1])
    
#                     mut_noise_dist = mutation_operator.mult_noise_distribution;
    
#                     mut_size_r = map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))
                    
#                     mut_size = abs.(reduce(vcat,mut_size_r))
                    
#                     if length(mut_size) != 0
#                         CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[3]) # new mult
#                     else
#                         mut_size = [0.]
#                     end
    
#                     hidedecorations!(ax1)
#                 else
#                     ax1 = Axis(plot_mut_hist[4,1])
    
#                     mut_noise_dist = mutation_operator.mult_noise_distribution;
    
#                     mut_size_r = map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,1,tr.H0-2),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))
                    
#                     mut_size = abs.(reduce(vcat,mut_size_r))
                    
#                     if length(mut_size) != 0
#                         CairoMakie.hist!(ax1,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[4]) # ex mult
#                     else
#                         mut_size = [0.]
#                     end
    
#                     hidedecorations!(ax1)
#                 end
    
#                 norm_pdf = [pdf(mut_noise_dist,t) for t in LinRange(min_t,max_t,100)];
    
#                 CairoMakie.lines!(ax1,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
#                 CairoMakie.vlines!(ax1,1.,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth) # ex all
    
#                 if type[1]
#                     ax2 = Axis(plot_mut_hist[3,2])
    
#                     mut_size_r = map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))
    
#                     mut_size = abs.(reduce(vcat,mut_size_r))
                    
#                     if length(mut_size) != 0
#                         CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[3]) # new mult
#                     else
#                         mut_size = [0.]
#                     end
                    
#                     hidedecorations!(ax2)
#                 else
#                     ax2 = Axis(plot_mut_hist[4,2])
    
#                     mut_size_r = map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))
    
#                     mut_size = abs.(reduce(vcat,mut_size_r))
                    
#                     if length(mut_size) != 0
#                         CairoMakie.hist!(ax2,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[4]) # ex mult
#                     else
#                         mut_size = [0.]
#                     end
    
#                     hidedecorations!(ax2)
#                 end

#                 norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];
    
#                 CairoMakie.lines!(ax2,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
#                 CairoMakie.vlines!(ax2,1.,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth) # ex all
    
#                 if type[1]
#                     ax3 = Axis(plot_mut_hist[3,3])
    
#                     mut_size_r = map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,tr.H0,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))
    
#                     mut_size = abs.(reduce(vcat,mut_size_r))
                    
#                     if length(mut_size) != 0
#                         CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[3]) # new mult
#                     else
#                         mut_size = [0.]
#                     end
    
#                     hidedecorations!(ax3) # ex all
#                 else
#                     ax3 = Axis(plot_mut_hist[4,3])
    
#                     mut_size_r = map(tr->get_mut_size_by_type_and_weight(tr,type,weight_id,tr.H0,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))
    
#                     mut_size = abs.(reduce(vcat,mut_size_r))
                    
#                     if length(mut_size) != 0
#                         CairoMakie.hist!(ax3,mut_size,bins = bins,normalization = :pdf,color = palette(:viridis, 4)[4]) # ex mult
#                     else
#                         mut_size = [0.]
#                     end
                    
#                     hidedecorations!(ax3) # ex all
#                 end
    
#                 norm_pdf = [pdf(mut_noise_dist,abs(t)) for t in LinRange(min_t,max_t,100)];
    
#                 CairoMakie.lines!(ax3,LinRange(min_t,max_t,100),norm_pdf,color = :red,linewidth = evo_config.wait_linewidth)
#                 CairoMakie.vlines!(ax3,1.,color = :red,linestyle = "--",linewidth = evo_config.wait_linewidth)
    
#                 push!(mult_ax,ax1)
#                 push!(mult_ax,ax2)
#                 push!(mult_ax,ax3)
#             end
    
#         end
    
#         linkxaxes!(additive_ax...)
#         linkxaxes!(mult_ax...)

#         colgap!(plot_mut_hist, 10)
#         rowgap!(plot_mut_hist, 10)

#         push!(weight_grid_layouts,plot_mut_hist)

#     end

#     for (label, layout) in zip(weight_names_latex, weight_grid_layouts[2:end])
#         Label(layout[1, 1, TopLeft()], label,
#             fontsize = 26,
#             font = :bold,
#             padding = (0, 5, 5, 0),
#             halign = :right)
#     end

#     Label(weight_grid_layouts[1][1, 1, TopLeft()], L"\text{All}",
#     fontsize = 26,
#     font = :bold,
#     padding = (0, 5, 5, 0),
#     halign = :right)

#     colors = palette(:viridis, 4)

# end

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

    if pred_type != :mss
        pop = filter(tr->tr.train_test_indicator == :train,trajectories_p)
    else
        pop = filter(tr->tr.train_test_indicator_mss == :train,trajectories_p)
    end

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


    if pred_type != :mss
        pop = filter(tr->tr.train_test_indicator == :test,trajectories_p)
    else
        pop = filter(tr->tr.train_test_indicator_mss == :test,trajectories_p)
    end

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

function create_prediction_summary!(fig,trajectories_p,pred_type,max_ce,vertex_to_predict_label,pred_config)

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

    if pred_type != :mss
        pop = filter(tr->tr.train_test_indicator == :train,trajectories_p)
    else
        pop = filter(tr->tr.train_test_indicator_mss == :train,trajectories_p)
    end

    mean_accuracies = []
    mean_null_accuracies = []
    roc = []

    for r in 1:max_ce
        
        pop_achieved = filter(tr->v_restricted_label_inclusion(tr,x->weight_edit_restriction_measure(x,r),:H0),pop)

        pop_not_achieved = filter(tr->!v_restricted_label_inclusion(tr,x->weight_edit_restriction_measure(x,r),:H0),pop)

        accuracies  = map(tr->v_restricted_accuracy(tr,x->weight_edit_restriction_measure(x,r),pred_type),pop_not_achieved)

        null_accuracies  = map(tr->v_restricted_accuracy(tr,x->weight_edit_restriction_measure(x,0),pred_type),pop_not_achieved)

        # if pred_type != :mss
        #     pop_na_labels = map(tr->vertex_to_predict_label[tr.inc_metagraph_vertices[tr.H0]],pop_not_achieved)
        #     pred_prob  = reduce(hcat,map(tr->v_restricted_probabilities(tr,x->weight_edit_restriction_measure(x,r),pred_type),pop_not_achieved)) |> transpose |> collect;
        #     roc_score = roc_auc_score(pop_na_labels,pred_prob,multi_class = "ovr", average = "weighted")
        #     push!(roc,roc_score)
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


    train_ac = CairoMakie.lines!(ax_train_accuracy,Float64.(mean_accuracies),color = :green, linestyle = "--",linewidth = pred_config.perf_linewidth)
    train_ac_null = CairoMakie.lines!(ax_train_accuracy,Float64.(mean_null_accuracies),color = :cyan,linestyle = "--",linewidth = pred_config.perf_linewidth)

    # if pred_type != :mss
    #     train_roc = CairoMakie.lines!(ax_train_accuracy,Float64.(roc),color = :purple,linestyle = "--",linewidth = pred_config.perf_linewidth)
    # end

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

    if pred_type != :mss
        pop = filter(tr->tr.train_test_indicator == :test,trajectories_p)
    else
        pop = filter(tr->tr.train_test_indicator_mss == :test,trajectories_p)
    end

    mean_accuracies = []
    mean_null_accuracies = []
    roc = []

    for r in 1:max_ce

        pop_achieved = filter(tr->v_restricted_label_inclusion(tr,x->weight_edit_restriction_measure(x,r),:H0),pop)

        pop_not_achieved = filter(tr->!v_restricted_label_inclusion(tr,x->weight_edit_restriction_measure(x,r),:H0),pop)

        accuracies  = map(tr->v_restricted_accuracy(tr,x->weight_edit_restriction_measure(x,r),pred_type),pop_not_achieved)

        null_accuracies  = map(tr->v_restricted_accuracy(tr,x->weight_edit_restriction_measure(x,0),pred_type),pop_not_achieved)

        # if pred_type != :mss
        #     pop_na_labels = map(tr->vertex_to_predict_label[tr.inc_metagraph_vertices[tr.H0]],pop_not_achieved)
        #     pred_prob  = reduce(hcat,map(tr->v_restricted_probabilities(tr,x->weight_edit_restriction_measure(x,r),pred_type),pop_not_achieved)) |> transpose |> collect
        #     roc_score = roc_auc_score(pop_na_labels,pred_prob,multi_class = "ovr", average = "weighted")
        #     push!(roc,roc_score)
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

    # if pred_type != :mss
    #     test_roc = CairoMakie.lines!(ax_test_accuracy,Float64.(roc),color = :purple,linestyle = "--",linewidth = pred_config.perf_linewidth)
    # end

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

function first_phenotypes!(fig,trajectories,n_weight,mut_prob,original_phenotype,all_first_pheno_n,thresh_afp,stripe_afp,left_st,evo_config)

    null_mut_n = [prob_k_mutations(k,mut_prob,10) for k in 1:n_weight];

    fm_num = map(tr->length(tr.mutant_info[1].mut_size),trajectories);
    fm_epi = map(tr->tr.epistasis[1],trajectories);
    s0_num = map(tr->length(tr.mutant_info[tr.H0-1].mut_size),trajectories);
    s0_epi = map(tr->tr.epistasis[tr.H0-1],trajectories);

    fig1 = fig[1,1] = GridLayout()
    fig2 = fig[1,2] = GridLayout()
    fig3 = fig[2,1] = GridLayout()
    fig4 = fig[2,2] = GridLayout()
    fig5 = fig[3,1] = GridLayout()
    fig6 = fig[3,2] = GridLayout()

    pheno_categories = [(.! thresh_afp) .& (.! stripe_afp) .& left_st, 
                        (.! thresh_afp) .& (.! stripe_afp) .& (.!left_st),
                        (thresh_afp) .& (.! stripe_afp) .& left_st,
                        (thresh_afp) .& (.! stripe_afp) .& (.!left_st),
                        stripe_afp,[true for _ in fm_num]]

    fig_subplots = [fig1,fig2,fig3,fig4,fig5,fig6]

    left_ind = [true,false,true,false,true,false]

    titles = [ L"\text{Phenotype: Left handed soft threshold }", 
                L"\text{Phenotype: Right handed soft threshold }",
                L"\text{Phenotype: Left handed hard threshold }", 
                L"\text{Phenotype: Right handed hard threshold }", 
                L"\text{Stripe}",L"\text{All}"]

    ylim = [0.,0.3,0.7]

    for (n,id) in enumerate(pheno_categories)

        subplot = fig_subplots[n]

        if n != length(pheno_categories)

            if n != length(pheno_categories)-1
                ax1 = Axis(subplot[1:2,:],title = titles[n],ygridvisible = false,xgridvisible = false)
                ax2 = Axis(subplot[3,:],ygridvisible = false,xgridvisible = false)
                ax3 = Axis(subplot[4,:],xlabel = L"\text{Weight changes per mutant}",ygridvisible = false,xgridvisible = false)

                for ph in all_first_pheno_n[id]
                    CairoMakie.lines!(ax1,ph, color = (:grey,0.1))
                end

                CairoMakie.lines!(ax1,original_phenotype, color = :red)

                epi_mutn_counts_lS0 = countmap(zip(fm_epi[id],fm_num[id]) |> collect)
                total_epi_mutn_lS0 = sum(values(epi_mutn_counts_lS0))
                epi_mutn_prop_lS0 = Dict(key=>value/total_epi_mutn_lS0 for (key,value) in epi_mutn_counts_lS0);
            
                grp_epi_lS0 = reduce(vcat,[[1,2,3,4] for i in 1:10])
                grp_epi_label_lS0 = reduce(vcat,[[:rse,:se,:ne,:sm] for i in 1:10])
                grp_mutn_lS0 = reduce(vcat,[[i,i,i,i] for i in 1:10])
            
                values_epi_mutn_lS0 = [haskey(epi_mutn_prop_lS0,key) ? epi_mutn_prop_lS0[key] : 0. for key in zip(grp_epi_label_lS0,grp_mutn_lS0)]

                color_epi_mutn_lS0 = [evo_config.pie_colors[i] for i in grp_epi_lS0];

                CairoMakie.barplot!(ax2,grp_mutn_lS0 ,values_epi_mutn_lS0, color = color_epi_mutn_lS0, stack = grp_epi_lS0)
                CairoMakie.lines!(ax2,null_mut_n, color = :black, linestyle = :dash )
                CairoMakie.scatter!(ax2,null_mut_n,color = :black,marker = 'x')

                ###################

                epi_mutn_counts_lS0 = countmap(zip(s0_epi[id],s0_num[id]) |> collect)
                total_epi_mutn_lS0 = sum(values(epi_mutn_counts_lS0))
                epi_mutn_prop_lS0 = Dict(key=>value/total_epi_mutn_lS0 for (key,value) in epi_mutn_counts_lS0);
            
                grp_epi_lS0 = reduce(vcat,[[1,2,3,4] for i in 1:10])
                grp_epi_label_lS0 = reduce(vcat,[[:rse,:se,:ne,:sm] for i in 1:10])
                grp_mutn_lS0 = reduce(vcat,[[i,i,i,i] for i in 1:10])
            
                values_epi_mutn_lS0 = [haskey(epi_mutn_prop_lS0,key) ? epi_mutn_prop_lS0[key] : 0. for key in zip(grp_epi_label_lS0,grp_mutn_lS0)]

                color_epi_mutn_lS0 = [evo_config.pie_colors[i] for i in grp_epi_lS0];

                CairoMakie.barplot!(ax3,grp_mutn_lS0,values_epi_mutn_lS0, color = color_epi_mutn_lS0, stack = grp_epi_lS0)
                CairoMakie.lines!(ax3,null_mut_n, color = :black, linestyle = :dash )
                CairoMakie.scatter!(ax3,null_mut_n,color = :black,marker = 'x')

                if left_ind[n]
                    CairoMakie.text!(ax1,75.,0.5, text = string(round(100*sum(id) / length(trajectories), digits = 2)) * "%")
                else
                    CairoMakie.text!(ax1,25.,0.3, text = string(round(100*sum(id) / length(trajectories), digits = 2)) * "%")
                end

                ax2.xticks = (1:n_weight,string.(1:n_weight))
                ax3.xticks = (1:n_weight,string.(1:n_weight))

                CairoMakie.xlims!(ax2,nothing,n_weight)
                CairoMakie.xlims!(ax3,nothing,n_weight)

                CairoMakie.ylims!(ax2,0,1.1*ylim[end])
                CairoMakie.ylims!(ax3,0,1.1*ylim[end])

                ax2.yticks = (ylim, string.(ylim))
                ax3.yticks = (ylim, string.(ylim))

                CairoMakie.text!(ax2,5,0.3, text = L"t=1")
                CairoMakie.text!(ax3,5,0.3, text = L"t=S_0")

                hidexdecorations!(ax2)
            else
                ax1 = Axis(subplot[1,:],title = titles[n],ygridvisible = false,xgridvisible = false)
                ax2 = Axis(subplot[2,:],xlabel = L"\text{Weight changes per mutant}",ygridvisible = false,xgridvisible = false)

                for ph in all_first_pheno_n[id]
                    CairoMakie.lines!(ax1,ph, color = (:grey,0.1))
                end

                CairoMakie.lines!(ax1,original_phenotype, color = :red)

                epi_mutn_counts_lS0 = countmap(zip(fm_epi[id],fm_num[id]) |> collect)
                total_epi_mutn_lS0 = sum(values(epi_mutn_counts_lS0))
                epi_mutn_prop_lS0 = Dict(key=>value/total_epi_mutn_lS0 for (key,value) in epi_mutn_counts_lS0);
            
                grp_epi_lS0 = reduce(vcat,[[1,2,3,4] for i in 1:10])
                grp_epi_label_lS0 = reduce(vcat,[[:rse,:se,:ne,:sm] for i in 1:10])
                grp_mutn_lS0 = reduce(vcat,[[i,i,i,i] for i in 1:10])
            
                values_epi_mutn_lS0 = [haskey(epi_mutn_prop_lS0,key) ? epi_mutn_prop_lS0[key] : 0. for key in zip(grp_epi_label_lS0,grp_mutn_lS0)]

                color_epi_mutn_lS0 = [evo_config.pie_colors[i] for i in grp_epi_lS0];

                CairoMakie.barplot!(ax2,grp_mutn_lS0 ,values_epi_mutn_lS0, color = color_epi_mutn_lS0, stack = grp_epi_lS0)
                CairoMakie.lines!(ax2,null_mut_n, color = :black, linestyle = :dash )
                CairoMakie.scatter!(ax2,null_mut_n,color = :black,marker = 'x')

                if left_ind[n]
                    CairoMakie.text!(ax1,75.,0.5, text = string(round(100*sum(id) / length(trajectories), digits = 2)) * "%")
                else
                    CairoMakie.text!(ax1,25.,0.3, text = string(round(100*sum(id) / length(trajectories), digits = 2)) * "%")
                end

                ax2.xticks = (1:n_weight,string.(1:n_weight))

                CairoMakie.xlims!(ax2,nothing,n_weight)
                CairoMakie.ylims!(ax2,0,1.1*ylim[end])

                ax2.yticks = (ylim, string.(ylim))

                CairoMakie.text!(ax2,5,0.3, text = L"t=1=S_0")
            end
        else
            ax1 = Axis(subplot[1,:],title = titles[n],ygridvisible = false,xgridvisible = false)
            ax2 = Axis(subplot[2,:],xlabel = L"\text{Weight changes per mutant}",ygridvisible = false,xgridvisible = false)

            epi_mutn_counts_lS0 = countmap(zip(fm_epi[id],fm_num[id]) |> collect)
            total_epi_mutn_lS0 = sum(values(epi_mutn_counts_lS0))
            epi_mutn_prop_lS0 = Dict(key=>value/total_epi_mutn_lS0 for (key,value) in epi_mutn_counts_lS0);
        
            grp_epi_lS0 = reduce(vcat,[[1,2,3,4] for i in 1:10])
            grp_epi_label_lS0 = reduce(vcat,[[:rse,:se,:ne,:sm] for i in 1:10])
            grp_mutn_lS0 = reduce(vcat,[[i,i,i,i] for i in 1:10])
        
            values_epi_mutn_lS0 = [haskey(epi_mutn_prop_lS0,key) ? epi_mutn_prop_lS0[key] : 0. for key in zip(grp_epi_label_lS0,grp_mutn_lS0)]

            color_epi_mutn_lS0 = [evo_config.pie_colors[i] for i in grp_epi_lS0];

            CairoMakie.barplot!(ax1,grp_mutn_lS0,values_epi_mutn_lS0, color = color_epi_mutn_lS0, stack = grp_epi_lS0)
            CairoMakie.lines!(ax1,null_mut_n, color = :black, linestyle = :dash )
            CairoMakie.scatter!(ax1,null_mut_n,color = :black,marker = 'x')

            ###################

            epi_mutn_counts_lS0 = countmap(zip(s0_epi[id],s0_num[id]) |> collect)
            total_epi_mutn_lS0 = sum(values(epi_mutn_counts_lS0))
            epi_mutn_prop_lS0 = Dict(key=>value/total_epi_mutn_lS0 for (key,value) in epi_mutn_counts_lS0);
        
            grp_epi_lS0 = reduce(vcat,[[1,2,3,4] for i in 1:10])
            grp_epi_label_lS0 = reduce(vcat,[[:rse,:se,:ne,:sm] for i in 1:10])
            grp_mutn_lS0 = reduce(vcat,[[i,i,i,i] for i in 1:10])
        
            values_epi_mutn_lS0 = [haskey(epi_mutn_prop_lS0,key) ? epi_mutn_prop_lS0[key] : 0. for key in zip(grp_epi_label_lS0,grp_mutn_lS0)]

            color_epi_mutn_lS0 = [evo_config.pie_colors[i] for i in grp_epi_lS0];

            CairoMakie.barplot!(ax2,grp_mutn_lS0,values_epi_mutn_lS0, color = color_epi_mutn_lS0, stack = grp_epi_lS0)
            CairoMakie.lines!(ax2,null_mut_n, color = :black, linestyle = :dash )
            CairoMakie.scatter!(ax2,null_mut_n,color = :black,marker = 'x')

            CairoMakie.xlims!(ax1,nothing,n_weight)
            CairoMakie.xlims!(ax2,nothing,n_weight)

            CairoMakie.ylims!(ax1,0,1.1*ylim[end])
            CairoMakie.ylims!(ax2,0,1.1*ylim[end])

            hidexdecorations!(ax1)

            ax1.xticks = (1:n_weight,string.(1:n_weight))
            ax2.xticks = (1:n_weight,string.(1:n_weight))

            ax1.yticks = (ylim, string.(ylim))
            ax2.yticks = (ylim, string.(ylim))

            CairoMakie.xlims!(ax1,0,n_weight)
            CairoMakie.xlims!(ax2,0,n_weight)

            CairoMakie.text!(ax1,5,0.3, text = L"t=1")
            CairoMakie.text!(ax2,5,0.3, text = L"t=S_0")

        end

        colgap!(subplot,Relative(0.05))
        rowgap!(subplot,Relative(0.05))
    end

    # colgap!(fig.layout,Relative(0.02))
    # rowgap!(fig.layout,Relative(0.01))

    labels_epi  = [L"\text{TD}",L"\text{SD}",L"\text{TI}",L"\text{SIC}"]

    symbol_epi = [PolyElement(color=c) for c in evo_config.pie_colors]

    Legend(fig[3, :, Bottom()], symbol_epi, labels_epi, framevisible=false,nbanks = 1,orientation = :horizontal,patchsize = (10, 10), colgap = 4, padding=(0.,0.,0.,100.))

end
