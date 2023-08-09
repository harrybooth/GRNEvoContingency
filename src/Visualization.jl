using Graphs
using GraphMakie
using DataInterpolations
using CairoMakie

function select_marker(edge_value)
    edge_value > 0 ? :ltriangle : :vline
end

function draw_grn(network,color_scheme)

    weight_indices = Tuple.(findall(ones(size(network)) .> 0));

    adjacency_matrix = vcat(network,zeros(1,4))

    ng = SimpleDiGraph(adjacency_matrix)

    edge_indices = [(src(i),dst(i)) for i in edges(ng)]

    edge_values = [round(adjacency_matrix[id...],digits = 3) for id in edge_indices]

    vertex_names = Dict(1=>"A",2=> "B", 3=> "C", 4=> "M")

    vertex_names = [vertex_names[i] for i in Graphs.vertices(ng)]

    fixed_layout(_) = [(0.5, 1.), (0., 0), (1.,0),(0.5,2.)]

    offsets = [[0.02, 0.0],
                [-0.03, -0.05],
                [0.01, -0.05],
                [0.0, 0.01]]

    e_colors = [color_scheme[findall(x->x==t,weight_indices)[1]] for t in edge_indices]

    f, ax, p = graphplot(ng,layout = fixed_layout,node_color = :black, elabels = string.(edge_values), nlabels = vertex_names,edge_color = e_colors,edge_width = 6.,arrow_size = 30.,arrow_shift = 0.1, arrow_attr = (; marker = [select_marker(edge_values[i]) for i in 1:ne(ng)]))

    p.elabels_rotation[] = Dict(i => edge_indices[i] == (1,4) ? 0.0 : Makie.automatic for i in 1:ne(ng))

    # offsets = 0.05 * (p[:node_pos][] .- p[:node_pos][][1])
    p.nlabels_offset[] = offsets

    autolimits!(ax)
    hidedecorations!(ax); hidespines!(ax)

    f,ax,p
end

function draw_grn!(ax,network,color_scheme)

    weight_indices = Tuple.(findall(ones(size(network)) .> 0));

    adjacency_matrix = vcat(network,zeros(1,4))

    ng = SimpleDiGraph(adjacency_matrix)

    edge_indices = [(src(i),dst(i)) for i in edges(ng)]

    edge_values = [round(adjacency_matrix[id...],digits = 3) for id in edge_indices]

    vertex_names = Dict(1=>"A",2=> "B", 3=> "C", 4=> "M")

    vertex_names = [vertex_names[i] for i in Graphs.vertices(ng)]

    fixed_layout(_) = [(0.5, 1.), (0., 0), (1.,0),(0.5,2.)]

    offsets = [[0.02, 0.0],
                [-0.075, -0.05],
                [0.025, -0.05],
                [0.0, 0.01]]

    e_colors = [color_scheme[findall(x->x==t,weight_indices)[1]] for t in edge_indices]

    graphplot!(ax,ng,layout = fixed_layout, node_color = :black, elabels = string.(edge_values), nlabels = vertex_names,edge_color = e_colors,edge_width = 6.,arrow_size = 30.,arrow_shift = 0.1, arrow_attr = (; marker = [select_marker(edge_values[i]) for i in 1:ne(ng)]), elabels_rotation = Dict(i => edge_indices[i] == (1,4) ? 0.0 : Makie.automatic for i in 1:ne(ng)),nlabels_offset = offsets)

    autolimits!(ax)
    hidedecorations!(ax); hidespines!(ax)
end

function draw_grn_cn!(ax,network)

    weight_indices = Tuple.(findall(ones(size(network)) .> 0));

    adjacency_matrix = vcat(network,zeros(1,4))

    ng = SimpleDiGraph(adjacency_matrix)

    edge_indices = [(src(i),dst(i)) for i in edges(ng)]

    edge_values = [round(adjacency_matrix[id...],digits = 3) for id in edge_indices]

    vertex_names = Dict(1=>"A",2=> "B", 3=> "C", 4=> "M")

    vertex_names = [vertex_names[i] for i in vertices(ng)]

    fixed_layout(_) = [(0.5, 1.), (0., 0), (1.,0),(0.5,2.)]

    offsets = [[0.03, 0.0],
                [-0.06, -0.06],
                [0.03, -0.06],
                [0.03, -0.06]]

    graphplot!(ax,ng,layout = fixed_layout, node_color = [:red,:blue,:orange,:green], node_size = 40., elabels = string.(edge_values), elabels_align = (:right,:center),nlabels = vertex_names,edge_color = :black,edge_width = 6.,arrow_size = 30.,arrow_shift = 0.1, arrow_attr = (; marker = [select_marker(edge_values[i]) for i in 1:ne(ng)]),selfedge_size = 0.25, nlabels_offset = offsets, elabels_rotation = Dict(i => edge_indices[i] == (1,4) ? 0.0 : Makie.automatic for i in 1:ne(ng)))

    autolimits!(ax)
    hidedecorations!(ax); hidespines!(ax)
end


function draw_grn_and_conc(w_ind,w_network,top_choice)

    fig = CairoMakie.Figure(resolution = (1800, 1000),fontsize = 22.)

    weight_indices = Tuple.(findall(ones(3,4) .> 0));

    vertex_names = Dict(1=>"A",2=> "B", 3=> "C", 4=> "M")

    line_colors = [:red,:blue,:orange,:green]

    ax1  = Axis(fig[1,1], backgroundcolor = "white", title = "Gene regulatory network : topology = " * top_choice )

    draw_grn_cn!(ax1,w_network)

    ax2  = Axis(fig[1,2], backgroundcolor = "white", title = "Half stripe",xlabel = "Tissue (cell id)", ylabel = "Concentration")

    for i in 1:3
        lines!(ax2,w_ind.phenotype.u[end][i,:],label = "Gene " * vertex_names[i], color = line_colors[i],linewidth = 7.)
    end

    lines!(ax2,map(x->morph(x),tissue),label = "Morphogen, M",linewidth = 7.,linestyle = "--",color = :green)

    leg = Legend(fig[1,3], ax2)

    return fig
end

function concat_matrices(v_m)

    r = zeros((size(v_m[1])...,length(v_m)))

    for i in 1:length(v_m)

        r[:,:,i] = v_m[i]

    end

    return r

end

function concat_matrices(v_m)

    r = zeros((size(v_m[1])...,length(v_m)))

    for i in 1:length(v_m)
        r[:,:,i] = v_m[i]
    end

    return r

end

function plot_clustered_geno_traj(geno_traj,fitness_traj,initial_fitness,initial_genotype,cluster_assignments,topology_name::String,n)

    # cluster_geno_trajectories!(gt,n_clusters,SqEuclidean());
    
    res = 1000

    n_clusters = length(unique(cluster_assignments))

    fig = CairoMakie.Figure(resolution = (res + 6.0*res/n_clusters, res),fontsize = 22.)

    color_scheme = cgrad(:tab20,categorical = true);

    weight_indices = Tuple.(findall(ones(3,4) .> 0));

    vertex_names = Dict(1=>"A",2=> "B", 3=> "C", 4=> "M")

    weight_names = [string(vertex_names[last(t)]) * "=>" * string(vertex_names[first(t)]) for t in weight_indices]

    gl = fig[1,2:4]

    grid_entries = Tuple.(findall(x-> x> 0,ones(Int64(ceil(sqrt(n_clusters))),Int64(ceil(sqrt(n_clusters))))))
    
    full_fitness_traj = collect(LinRange(initial_fitness,1.,n))

    ax_list = []

    full_geno_traj = []

    for i in 1:length(fitness_traj)
        full_weight_traj = zeros(size(geno_traj[i],1),length(LinRange(initial_fitness,1.,n)))
        for wi in 1:12
            unique_geno_visits = unique(hcat(fitness_traj[i],geno_traj[i][wi,:]), dims = 1)
            fitness_timestamps = unique_geno_visits[:,1]
            itp_g = DataInterpolations.ConstantInterpolation(unique_geno_visits[:,2],fitness_timestamps);
            full_weight_traj[wi,:] = [itp_g(t) for t in LinRange(initial_fitness,1.,n)]
        end
        push!(full_geno_traj,full_weight_traj)
    end

    # fgt_data = hcat(map(x->vec(x),full_geno_traj)...)

    for (n,cl) in enumerate(unique(cluster_assignments))

        ass_cl = full_geno_traj[findall(x->x==cl,cluster_assignments)];

        all_m = concat_matrices(ass_cl)

        mt = mean(all_m,dims = 3)[:,:,1]
        std_t = std(all_m,dims = 3)[:,:,1]

        ax1  = Axis(gl[grid_entries[n]...], backgroundcolor = "white", xlabel = "Fitness", ylabel = "Weight value", title = "Final mechanism: " * string(cl))

        for i in 1:12
            CairoMakie.lines!(ax1,full_fitness_traj,mt[i,:], label = "Weight " * weight_names[i],color = color_scheme[i],linewidth = 5.)
            # CairoMakie.band!(ax1,full_fitness_traj,mt[i,:] .- std_t[i,:],mt[i,:] .+ std_t[i,:],color = (color_scheme[i],0.5))
        end

        push!(ax_list,ax1)

    end

    cm = countmap(cluster_assignments)

    cd = zeros(n_clusters)

    for (n,cl) in enumerate(unique(cluster_assignments))
        cd[n] = cm[cl] / length(cluster_assignments)
    end
    
    gr = fig[1,1]

    leg = Legend(gr, ax_list[1])

    # Barplot

    gbp = fig[1,5]

    ocd_labels = unique(cluster_assignments)[sortperm(cd)]

    # ocd = cd[ocd_labels]

    ax2  = Axis(gbp[1,1], backgroundcolor = "white", xlabel = "Mechanism", ylabel = "Frequency",xticks = (1:length(ocd_labels),string.(ocd_labels)),title = "Mechanism association frequency",xticklabelrotation=45)

    CairoMakie.barplot!(ax2,sort(cd))

    # GRN

    ax3  = Axis(gbp[2,1], backgroundcolor = "white",title = topology_name * " - start network")

    draw_grn!(ax3,reshape(initial_genotype,(3,4)),color_scheme)

    return fig
end

# function plot_clustered_geno_traj(gt::GenoTrajectories, cluster_assignments, topology_name::String,n)

#     # cluster_geno_trajectories!(gt,n_clusters,SqEuclidean());
    
#     res = 1000

#     n_clusters = length(unique(cluster_assignments))

#     fig = CairoMakie.Figure(resolution = (res + 6.0*res/n_clusters, res),fontsize = 22.)

#     color_scheme = cgrad(:tab20,categorical = true);

#     weight_indices = Tuple.(findall(ones(3,4) .> 0));

#     vertex_names = Dict(1=>"A",2=> "B", 3=> "C", 4=> "M")

#     weight_names = [string(vertex_names[last(t)]) * "=>" * string(vertex_names[first(t)]) for t in weight_indices]

#     gl = fig[1,2:4]

#     grid_entries = Tuple.(findall(x-> x> 0,ones(Int64(ceil(sqrt(n_clusters))),Int64(ceil(sqrt(n_clusters))))))
    
#     full_fitness_traj = collect(LinRange(gt_mi.initial_fitness,1.,n))

#     ax_list = []

#     full_geno_traj = []

#     for i in 1:length(gt.fitness_traj)
#         full_weight_traj = zeros(size(gt.geno_traj[i],1),length(LinRange(gt_mi.initial_fitness,1.,n)))
#         for wi in 1:12
#             unique_geno_visits = unique(hcat(gt.fitness_traj[i],gt.geno_traj[i][wi,:]), dims = 1)
#             fitness_timestamps = unique_geno_visits[:,1]
#             itp_g = DataInterpolations.ConstantInterpolation(unique_geno_visits[:,2],fitness_timestamps);
#             full_weight_traj[wi,:] = [itp_g(t) for t in LinRange(gt_mi.initial_fitness,1.,n)]
#         end
#         push!(full_geno_traj,full_weight_traj)
#     end

#     # fgt_data = hcat(map(x->vec(x),full_geno_traj)...)

#     for cl in unique(cluster_assignments)

#         ass_cl = full_geno_traj[findall(x->x==cl,cluster_assignments)];

#         all_m = concat_matrices(ass_cl)

#         mt = mean(all_m,dims = 3)[:,:,1]
#         std_t = std(all_m,dims = 3)[:,:,1]

#         ax1  = Axis(gl[grid_entries[cl]...], backgroundcolor = "white", xlabel = "Fitness", ylabel = "Weight value", title = "Cluster " * string(cl))

#         for i in 1:12
#             CairoMakie.lines!(ax1,full_fitness_traj,mt[i,:], label = "Weight " * weight_names[i],color = color_scheme[i],linewidth = 5.)
#         #     CairoMakie.band!(ax1,full_fitness_traj,mt[i,:] .- std_t[i,:],mt[i,:] .+ std_t[i,:],color = (color_scheme[i],0.5))
#         end

#         push!(ax_list,ax1)

#     end

#     cm = countmap(cluster_assignments)

#     cd = zeros(n_clusters)

#     for cl in unique(cluster_assignments)
#         cd[cl] = cm[cl] / length(cluster_assignments)
#     end
    
#     gr = fig[1,1]

#     leg = Legend(gr, ax_list[1])

#     # Barplot

#     gbp = fig[1,5]

#     ocd_labels = sortperm(cd)

#     ocd = cd[ocd_labels]

#     ax2  = Axis(gbp[1,1], backgroundcolor = "white", xlabel = "Cluster", ylabel = "Frequency",xticks = (1:length(ocd_labels),string.(ocd_labels)),title = "Cluster association frequency")

#     CairoMakie.barplot!(ax2,ocd)

#     # GRN

#     ax3  = Axis(gbp[2,1], backgroundcolor = "white",title = topology_name * " - start network")

#     draw_grn!(ax3,reshape(gt.initial_genotype,(3,4)),color_scheme)

#     return fig
# end

function draw_fitness_slices(ll)
    fig = CairoMakie.Figure(resolution = (2000,1800),fontsize = 30.)

    color_scheme = palette(:tab20)

    weight_indices = Tuple.(findall(viable_mutations .> 0));

    vertex_names = Dict(1=>"A",2=> "B", 3=> "C", 4=> "M")

    weight_names = [string(vertex_names[last(t)]) * "=>" * string(vertex_names[first(t)]) for t in weight_indices]
        
    ax_list = []

    for wi in 1:length(weight_indices)

        ax1  = Axis(fig[weight_indices[wi]...], backgroundcolor = "white", xlabel="Mutation size", ylabel="Fitness", title = weight_names[wi]) 

        # CairoMakie.scatter!(ax1,phate_pheno_all,markersize = 3.5,color = map(x->x == wi ? color_scheme[wi] : :grey,traj_fm_class_v[choice_full_v]))

        CairoMakie.lines!(ax1,ll.sample_points,ll.slice_fitnesses[weight_indices[wi]...,:],color = color_scheme[wi])

        CairoMakie.scatter!(ax1,[(0.,ll.origin_fitness)],color = color_scheme[wi],markersize = 15.)

        push!(ax_list,ax1)

    end

    linkxaxes!(ax_list...)
    linkyaxes!(ax_list...)

    # CairoMakie.lines(ll.sample_points,ll.slice_fitnesses[(1,3)...,:])

    return fig
end

function draw_fitness_slices_scan(ll)
    fig = CairoMakie.Figure(resolution = (2000,1800),fontsize = 30.)

    color_scheme = palette(:tab20)

    weight_indices = Tuple.(findall(viable_mutations .> 0));

    vertex_names = Dict(1=>"A",2=> "B", 3=> "C", 4=> "M")

    weight_names = [string(vertex_names[last(t)]) * "=>" * string(vertex_names[first(t)]) for t in weight_indices]
        
    ax_list = []

    for wi in 1:length(weight_indices)

        ax1  = Axis(fig[weight_indices[wi]...], backgroundcolor = "white", xlabel="Weight Value", ylabel="Fitness", title = weight_names[wi]) 

        # CairoMakie.scatter!(ax1,phate_pheno_all,markersize = 3.5,color = map(x->x == wi ? color_scheme[wi] : :grey,traj_fm_class_v[choice_full_v]))

        CairoMakie.lines!(ax1,ll.sample_points[weight_indices[wi]...,:],ll.slice_fitnesses[weight_indices[wi]...,:],color = color_scheme[wi])

        CairoMakie.scatter!(ax1,[(Float64(ll.origin.genotype.p[1][weight_indices[wi]...]),ll.origin_fitness)],color = color_scheme[wi],markersize = 15.)

        push!(ax_list,ax1)

    end

    linkxaxes!(ax_list...)
    linkyaxes!(ax_list...)

    # CairoMakie.lines(ll.sample_points,ll.slice_fitnesses[(1,3)...,:])

    return fig
end

function create_bar_from_tuples!(ax,data,color_dict,numeric_x = true)
    cat = unique(first.(data))

    if numeric_x
        all_values = sort(unique(last.(data)))
    else
        all_values = unique(last.(data))
    end

    chart_data = []
    color_data = []

    for (n,c) in enumerate(cat)
        n_cat = length(filter(x->first(x) == c,data))
        for v in all_values
            n_val = length(filter(x->(first(x) == c) & (last(x) == v) ,data))
            push!(chart_data, [n,v,n_val/n_cat])
            push!(color_data, color_dict[c])
        end
    end

    cd = reduce(hcat,chart_data) 

    CairoMakie.barplot!(ax,cd[2,:],cd[3,:], dodge = Int.(cd[1,:]),color = color_data)

    if numeric_x
        ax.xticks = (1:length(all_values),string.(all_values))
    end

end

function pairplot_makie_ass(X,ass,colorscheme)

    fig = CairoMakie.Figure(resolution = (5000,5000),fontsize = 40.)

    vertex_names = Dict(1=>"A",2=> "B", 3=> "C", 4=> "M")

    weight_indices = Tuple.(findall(ones(3,4) .> 0));

    weight_names = [string(vertex_names[last(t)]) * "=>" * string(vertex_names[first(t)]) for t in weight_indices]

    color_scheme = palette(:tab10)

    n_weight = 10

    grid_entries = Tuple.(findall(x->x > 0, ones(n_weight,n_weight)))

    ax_list = []

    for entry in grid_entries

        # if (entry[1] <= 12) && (entry[2] <= 12)
        if entry[1] == entry[2]
            ax = Axis(fig[entry...], xlabel = weight_names[entry[2]],ylabel = weight_names[entry[1]])
            CairoMakie.density!(ax,X[entry[2],:])

            if entry[2] != 1
                hideydecorations!(ax)
            end
        
            if (entry[1] != n_weight)
                hidexdecorations!(ax)
            end

            push!(ax_list,ax)

        elseif entry[1] >= entry[2]
            ax = Axis(fig[entry...], xlabel = weight_names[entry[2]],ylabel = weight_names[entry[1]])
            CairoMakie.scatter!(ax,X[[entry[2],entry[1]],:],color = ass,colormap  = colorscheme)

            if entry[2] != 1
                hideydecorations!(ax)
            end
        
            if (entry[1] != n_weight)
                hidexdecorations!(ax)
            end

            push!(ax_list,ax)
            
        else
            nothing
        end

    end

    rowgap!(fig.layout, 2.)
    colgap!(fig.layout, 4.)

    linkyaxes!(ax_list...)
    linkxaxes!(ax_list...)

    return fig

end