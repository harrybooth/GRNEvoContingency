function plot_convergence_rate!(ax,conv,conv_time)

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
    
    draw_grn!(ax_geno,orig_network,draw_config,node_colors,fs,false)

    hidedecorations!(ax)

end

function plot_dynamical_summary!(fig,trajectories,embedding,top_n,minimal_motif_count,sorted_uep,sorted_counts_uep,uep_position_dict,vertex_top_map,draw_config,example_mst,tr_choice,fontsize)

    trajectories_p_d = filter(tr->tr.inc_metagraph_vertices[end] ∈ sorted_uep[1:top_n],trajectories);

    mo_umap = fig[1:4, 1:3] = GridLayout()

    ex1 = fig[1:3, 4:7] = GridLayout()
    rmh0 = fig[4, 4:7] = GridLayout()

    top_n_dict = Dict(v_id=>pos for (pos,v_id) in enumerate(sorted_uep[1:top_n]))

    #### Motif Distribution

    color_sorted_counts_uep = [i <= top_n ? color_scheme[i] : :grey for i in 1:length(sorted_counts_uep)]

    view_sorted_uep_id = sorted_counts_uep .> minimal_motif_count

    other = sum(sorted_counts_uep[.!view_sorted_uep_id])

    n_other = sum(.!view_sorted_uep_id)

    n_norm = sum(view_sorted_uep_id)

    view_sorted_uep_counts = vcat(sorted_counts_uep[view_sorted_uep_id],[other])

    view_color_sorted_counts_uep = vcat(color_sorted_counts_uep[view_sorted_uep_id],[:grey])

    axis_label = vcat(collect(string.(1:length(view_sorted_uep_counts)-1)),["other"])

    ##############

    ax1 = Axis(mo_umap[1:2,1:top_n],title = "Top " * string(top_n) * " MST : " * string(sum(sorted_counts_uep[1:top_n])) * " trajectories", xlabel = L"\text{Dynamics: UMAP 1}", ylabel = L"\text{Dynamics: UMAP 2}")

    CairoMakie.scatter!(ax1,embedding, color = [haskey(top_n_dict,i) ? (color_scheme[top_n_dict[i]],0.5) : (:grey,0.5) for i in end_parents],markersize = 8.)

    hidedecorations!(ax1, label = false)

    for i in 1:top_n

        ax_geno = Axis(mo_umap[3,i], backgroundcolor = (color_scheme[i],0.5))

        top = Int.(reshape(vertex_top_map[sorted_uep[i]],(3,4)))

        # draw_grn_layout!(ax_geno,top,e_width,vertex_size,arrow_size,arrow_shift,sw,fixed_layout,selfedge_size,node_colors,false)

        draw_grn!(ax_geno,vertex_top_map[sorted_uep[i]],draw_config,node_colors,fontsize,false)
    end

    #######################

    ax_mo = Axis(mo_umap[4,1:top_n],ylabel  = L"\text{% of trajectories}")

    CairoMakie.barplot!(ax_mo,view_sorted_uep_counts ./ sum(view_sorted_uep_counts),color = view_color_sorted_counts_uep)

    # CairoMakie.text!(ax,[Point2f(n_norm+0.4,1.2*other/sum(view_sorted_uep_counts))],text = ["# MN-other = " * string(n_other)],fontsize = fontsize)

    ax_mo.xticks = (1:length(view_sorted_uep_counts),"MST " .* axis_label)

    n_other = sum(.!view_sorted_uep_id)

    ax_mo.xticklabelrotation = 45.

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

    n_unique_top = length(tr_traj)

    # tr_fd = data["fitness_traj"][conv][tr_data_id][tr_choice];

    tr_fd = create_full_fitness_traj(tr_data[tr_choice].fitness_traj_tuple,tr_data[tr_choice].wait_times)

    tr_top_fitness_id = uniqueid(tr_fd)[tr_traj_id]

    tr_fd_coarse = map(x->x[1],tr_fd)
    tr_fd_refine = map(x->x[2],tr_fd)

    tr_top_fitness_sc = [(x,tr_fd_refine[x]) for x in tr_top_fitness_id]

    ###################

    ax_fitness = Axis(ex1[1:2,1:length(tr_phenotypes)],xlabel = L"\text{Generation}")

    ax_fitness_c = Axis(ex1[1:2,1:length(tr_phenotypes)], yticklabelcolor = :blue, yaxisposition = :right, ylabel = L"\mathcal{F}_{C}(g_c)", ylabelcolor = :blue, ylabelrotation = 2pi)

    hidespines!(ax_fitness_c)
    hideydecorations!(ax_fitness_c,label = false,ticklabels = false,ticks = false,minorticks = false)
    hidexdecorations!(ax_fitness_c)

    hideydecorations!(ax_fitness,label = false,ticklabels = false,ticks = false,minorticks = false)

    CairoMakie.lines!(ax_fitness,tr_fd_refine, color = :grey, linewidth = 3.)
    CairoMakie.lines!(ax_fitness_c,tr_fd_coarse, linestyle = "--", color = :blue,linewidth = 3.)

    CairoMakie.scatter!(ax_fitness,tr_top_fitness_sc, color = [palette(:tab10)[i+1] for i in 1:length(tr_top_fitness_sc)], markersize = 20.)

    ####################

    for i in 1:length(tr_phenotypes)

        if i >=3
            ax_geno = Axis(ex1[3:4,i], backgroundcolor = (color_scheme[example_mst],0.5))
        else
            ax_geno = Axis(ex1[3:4,i], backgroundcolor = RGBf(0.98, 0.98, 0.98))
        end

        ax_pheno = Axis(ex1[5,i])

        for g in 1:3
            CairoMakie.lines!(ax_pheno,tr_phenotypes[i][g,:],linewidth = 4., color = color_scheme[4:6][g])
        end

        CairoMakie.hidedecorations!(ax_pheno)

        # top = reshape(tr_traj[i],(3,4))

        # draw_grn_layout!(ax_geno,top,e_width,vertex_size,arrow_size,arrow_shift,sw,fixed_layout,selfedge_size,node_colors,false)

        draw_grn!(ax_geno,tr_traj[i],draw_config,node_colors,fontsize,false)

    end

    ##########################

    all_prop  = []
    all_dodge  = []
    all_x  = []

    for n in 1:top_n

        pop = filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories_p_d)

        pop_equal = filter(tr->tr.minimal_stripe_subgraphs[tr.H0] == tr.minimal_stripe_subgraphs[end], pop)

        pop_incl = filter(tr->Bool(test_inclusion(tr.minimal_stripe_subgraphs[end],tr.minimal_stripe_subgraphs[tr.H0])),pop)

        n_pop = length(pop)

        proportions = [length(pop_equal)/n_pop,(length(pop_incl) - length(pop_equal))/n_pop,(length(pop) - length(pop_incl))/n_pop]

        x = [1,2,3]

        dodge = [n,n,n]

        push!(all_prop,proportions)
        push!(all_dodge,dodge)
        push!(all_x,x)

    end

    ax_rh0 = Axis(rmh0[1,1])

    x = reduce(vcat,all_x)
    dodge = reduce(vcat,all_dodge)
    proportions = reduce(vcat,all_prop)

    CairoMakie.barplot!(ax_rh0,x,proportions,color = [color_scheme[n] for n in dodge],dodge = dodge)

    CairoMakie.hidedecorations!(ax_rh0,label = false,ticklabels = false,ticks = false,minorticks = false)

    ax_rh0.xticks = (1:3,[L"M^{(i)}_{H_{0}} = M^{(i)}_{N_i}",L"M^{(i)}_{H_{0}} \subset M^{(i)}_{N_i}",L"M^{(i)}_{H_{0}} \nsubset M^{(i)}_{N_i}"])

    for (label, layout) in zip(["A", "B", "C"], [mo_umap, ex1, rmh0])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = 22,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right)
    end

end

function create_evo_summary!(fig,trajectories,top_n,mutation_operator,sorted_uep, vertex_top_map, draw_config, node_colors,fontsize,color_scheme)

    all_wait_times = reduce(hcat,[cumulative_wait_time(tr) for tr in trajectories]);

    colors = reverse(palette(:tab10)[1:4])

    ax_wait_list = []

    scale = 10^4

    name = L"(10^4)"

    ax_wait_list = []
    ax_wait_2_list = []

    wt_l_list = []
    wt_s_list = []

    min_t_n = -3.
    max_t_n = 3.

    min_t_u = -mutation_operator.max_w
    max_t_u = mutation_operator.max_w

    for n in 1:top_n

        plot_geno = fig[n, 1] = GridLayout()
        plot_wait = fig[n, 2] = GridLayout()
        plot_mut_hist = fig[n, 3:4] = GridLayout()
        plot_epi_types = fig[n, 5:6] = GridLayout()

        if n==1
            ax_geno = Axis(plot_geno[1,1],backgroundcolor = (color_scheme[n],0.5),title =L"\text{Minimal Stripe Topology}")
        else
            ax_geno = Axis(plot_geno[1,1],backgroundcolor = (color_scheme[n],0.5))
        end

        # top = Int.(reshape(vertex_top_map[sorted_uep[n]],(3,4)))

        # draw_grn_layout!(ax_geno,top,e_width,vertex_size,arrow_size,arrow_shift,sw,fixed_layout,selfedge_size,node_colors,false)

        draw_grn!(ax_geno,vertex_top_map[sorted_uep[n]],draw_config,node_colors,fontsize,false)

        if n == 1
            ax_wait = Axis(plot_wait[1,1], title = L"\mathbb{E}[\text{total weight edits}]")
        else
            ax_wait = Axis(plot_wait[1,1])
        end

        ax_wait_2 = Axis(plot_wait[1,1], yticklabelcolor = :red, yaxisposition = :right,yscale = log10)

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

        ax_wait.xticks = (1:3,[L"t<H_{0}",L"t=H_{0}",L"t>H_{0}" ])
        
        #format y ticks to latex numbers

        CairoMakie.hidexdecorations!(ax_wait,label = false,ticklabels = false,ticks = false,minorticks = false)
        CairoMakie.hideydecorations!(ax_wait,label = false,ticklabels = false,ticks = false,minorticks = false,grid = false)

        push!(ax_wait_list,ax_wait)

        ############################

        mean_wait = mean(all_wait_times[:,findall(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)],dims = 2)[:,1]
        mean_wait_v = reduce(vcat,mean_wait)

        # mean_wait_v = log.(reduce(vcat,mean_wait))

        mean_wait_type_labels = [1,2,3]

        wt_l = CairoMakie.lines!(ax_wait_2,mean_wait_type_labels,mean_wait_v,color = :red)
        wt_s = CairoMakie.scatter!(ax_wait_2,mean_wait_type_labels,mean_wait_v,color = :red)

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

        CairoMakie.pie!(ax_epi_lH0,epi_counts_prop,radius = 4,color = colors,
        inner_radius = 2,
        strokecolor = :white,
        strokewidth = 5)

        CairoMakie.hidedecorations!(ax_epi_lH0)

        if n == 1
            ax_epi_H0 = Axis(plot_epi_types[1,2],title = L"t=H_{0}")
        else
            ax_epi_H0 = Axis(plot_epi_types[1,2])
        end

        epi_counts = map(tr->tr.other[tr.H0-1],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_H0,epi_counts_prop,radius = 4,color = colors,
        inner_radius = 2,
        strokecolor = :white,
        strokewidth = 5)

        CairoMakie.hidedecorations!(ax_epi_H0)

        if n == 1
            ax_epi_uH0 = Axis(plot_epi_types[1,3],title = L"t>H_{0}")
        else
            ax_epi_uH0 = Axis(plot_epi_types[1,3])
        end

        epi_counts = reduce(vcat,map(tr->tr.other[tr.H0:end],filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

        epi_counts_prop = calculate_epi_class_proportion(epi_counts)

        CairoMakie.pie!(ax_epi_uH0,epi_counts_prop,radius = 4,color = colors,
        inner_radius = 2,
        strokecolor = :white,
        strokewidth = 5)

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

            pv1 = pvalue(OneSampleADTest(Float64.(mut_size),mut_noise_dist))

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

            CairoMakie.lines!(ax1,LinRange(min_t,max_t,100),norm_pdf,color = :red)
            CairoMakie.vlines!(ax1,0,color = :red,linestyle = "--")

            mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,tr.H0-1,tr.H0-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

            pv2 = pvalue(OneSampleADTest(Float64.(mut_size),mut_noise_dist))

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

            CairoMakie.lines!(ax2,LinRange(min_t,max_t,100),norm_pdf,color = :red)
            CairoMakie.vlines!(ax2,0,color = :red,linestyle = "--")

            mut_size = reduce(vcat,map(tr->get_mut_size_by_type(tr,type,tr.H0+1,length(tr.geno_traj)-1),filter(tr->tr.inc_metagraph_vertices[end] == sorted_uep[n],trajectories)))

            pv3 = pvalue(OneSampleADTest(Float64.(mut_size),mut_noise_dist))

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
        colgap!(plot_epi_types, 10)
        rowgap!(plot_epi_types, 10)
    end

    labels = [L"\text{RSE}",L"\text{Sign epistasis}",L"\text{No epistasis}",L"\text{Single mutation}"]

    Legend(fig[top_n+1, 5:6], [PolyElement(color=c) for c in colors], labels, framevisible=false,nbanks = 2,orientation = :horizontal)

    colors = palette(:viridis, 3)

    Legend(fig[top_n+1,2:4],
        vcat([[wt_s_list[1], wt_l_list[1]]],[PolyElement(color=c) for c in colors]),
        vcat([L"\mathbb{E}[\text{time to accept}]"],[L"\text{weight edits : existing}",L"\text{weight edits : new}",L"\text{weight edits : remove}"]),framevisible=false,nbanks = 2,orientation = :horizontal)


    linkyaxes!(ax_wait_list...)
    linkyaxes!(ax_wait_2_list...)

end

function create_weight_edit_summary!(fig,n,trajectories,mutation_op,sorted_uep, vertex_top_map, draw_config, node_colors,fontsize,color_scheme)

    grid_values = Tuple.(findall(ones(3,4) .> 0))

    colors = reverse(palette(:tab10)[1:4])

    ax_wait_list = []

    scale = 10^4

    name = L"(10^4)"

    ax_wait_list = []
    ax_wait_2_list = []

    wt_l_list = []
    wt_s_list = []

    min_t_n = -3.
    max_t_n = 3.

    min_t_u = -max_w
    max_t_u = max_w

    plot_geno = fig[grid_values[11]...] = GridLayout()

    ax_geno = Axis(plot_geno[1,1],backgroundcolor = (color_scheme[n],0.5),title =L"\text{Minimal Stripe Topology}")

    # top = Int.(reshape(vertex_top_map[sorted_uep[n]],(3,4)))

    # draw_grn_layout!(ax_geno,top,e_width,vertex_size,arrow_size,arrow_shift,sw,fixed_layout,selfedge_size,node_colors,true)

    draw_grn!(ax_geno,vertex_top_map[sorted_uep[n]],draw_config,node_colors,fontsize,false)

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
            ax_geno = Axis(plot_geno[1,1],backgroundcolor = (color_scheme[n],0.5),title =L"\text{Minimal Stripe Topology}")
        else
            ax_geno = Axis(plot_geno[1,1],backgroundcolor = (color_scheme[n],0.5))
        end

        # top = Int.(reshape(vertex_top_map[sorted_uep[n]],(3,4)))

        # draw_grn_layout!(ax_geno,top,e_width,vertex_size,arrow_size,arrow_shift,sw,fixed_layout,selfedge_size,node_colors,false)

        draw_grn!(ax_geno,vertex_top_map[sorted_uep[n]],draw_config,node_colors,fontsize,false)

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