mutable struct DrawGRNConfig

    l 

    shift_factor_node
    shift_factor_edge
    shift_factor_label

    arrow_size_act
    arrow_size_inh
    arrow_attr

    self_arc_radius

    node_size
    A_pos
    B_pos
    C_pos
    M_pos
    fixed_layout_scatter

    shift_axis_edge 
    shift_axis_edge_m
end

function DrawGRNConfig(l,shift_factor_node,shift_factor_edge,shift_factor_label,arrow_size_act,arrow_size_inh,arrow_attr,self_arc_radius,node_size,A_pos,B_pos,C_pos,M_pos,shift_axis_edge,shift_axis_edge_m)

    DrawGRNConfig(l,shift_factor_node,shift_factor_edge,shift_factor_label,arrow_size_act,arrow_size_inh,arrow_attr,self_arc_radius,node_size,A_pos,B_pos,C_pos,M_pos,[A_pos, B_pos, C_pos,M_pos],shift_axis_edge,shift_axis_edge_m)

end

function fs18_default()

    fs = 18.

    l = fs / 10

    shift_factor_node =  0.85

    shift_factor_edge(x) = x ? 0.05*l : 0.

    shift_factor_label = 0.15*shift_factor_node

    arrow_size_act = 0.8*fs
    arrow_size_inh = 1.5*fs

    arrow_attr = (;linewidth = 3., color = :black)

    self_arc_radius = l/7

    node_size = 1.3*fs

    B_pos = Point2f((0.,0.))
    A_pos = Point2f((-l*cos(pi/3),l*sin(pi/3)))
    C_pos = Point2f((l*cos(pi/3),l*sin(pi/3)))

    M_pos = Point2f((A_pos[1],A_pos[2] + l/2))

    DrawGRNConfig(l,shift_factor_node,shift_factor_edge,shift_factor_label,arrow_size_act,arrow_size_inh,arrow_attr,self_arc_radius,node_size,A_pos,B_pos,C_pos,M_pos,[A_pos, B_pos, C_pos,M_pos],l/3,l/5)

end

function fs12_default()

    fs = 12.

    l = fs / 10

    # shift_factor_node =  l/(2.5)

    shift_factor_node = 0.8

    shift_factor_edge(x) = x ? 0.075*l : 0.

    shift_factor_label = 0.15*shift_factor_node

    arrow_size_act = 0.4*fs
    arrow_size_inh = 0.8*fs

    arrow_attr = (;linewidth = 1., color = :black)

    self_arc_radius = l/6

    node_size = 0.7*fs

    B_pos = Point2f((0.,0.))
    A_pos = Point2f((-l*cos(pi/3),l*sin(pi/3)))
    C_pos = Point2f((l*cos(pi/3),l*sin(pi/3)))

    M_pos = Point2f((A_pos[1],A_pos[2] + l/2))

    DrawGRNConfig(l,shift_factor_node,shift_factor_edge,shift_factor_label,arrow_size_act,arrow_size_inh,arrow_attr,self_arc_radius,node_size,A_pos,B_pos,C_pos,M_pos,[A_pos, B_pos, C_pos,M_pos],l/2,l/3)

end

function plot_shifted_arrows!(ax,arrow_type,pos_start::Point2f,pos_end::Point2f,shift_factor)
    dir_travel = pos_end - pos_start

    # n_dir_travel = dir_travel / norm(dir_travel)

    pos_start_new = pos_start + shift_factor*dir_travel
    pos_end_new = pos_start + (1-shift_factor)*dir_travel

    dir_travel_new = pos_end_new - pos_start_new

    CairoMakie.arrows!(ax,[pos_start_new[1]],[pos_start_new[2]],[dir_travel_new[1]],[dir_travel_new[2]], arrowhead = arrow_type)
end

function plot_shifted_arrows!(ax,arrow_type,arrow_size,pos_start::Point2f,pos_end::Point2f,annotation,shift_factor_node, shift_factor_edge,shift_factor_label,arrow_attr)
    dir_travel = pos_end - pos_start

    norm_dt = Point2f((dir_travel[2],-dir_travel[1])) / norm(dir_travel)

    # n_dir_travel = dir_travel / norm(dir_travel)

    pos_start_new = (pos_start + (1-shift_factor_node)*dir_travel) + shift_factor_edge*norm_dt
    pos_end_new = (pos_start + shift_factor_node*dir_travel) + shift_factor_edge*norm_dt

    dir_travel_new = pos_end_new - pos_start_new

    mid_point = pos_start_new + 0.5*dir_travel_new

    mid_point_text = mid_point + shift_factor_label*norm(dir_travel)*norm_dt

    CairoMakie.arrows!(ax,[pos_start_new[1]],[pos_start_new[2]],[dir_travel_new[1]],[dir_travel_new[2]]; arrowhead = arrow_type,arrowsize = arrow_size,arrow_attr...)

    if ! isnothing(annotation)
        CairoMakie.text!(ax,mid_point_text,text = annotation, align = (:center,:baseline),color = arrow_attr.color)
    end
end

function plot_self_arc!(ax,arrow_type,arrow_size,radius,pos_start::Point2f,annotation,node_id,arrow_attr)

    if node_id == :B

        centre = pos_start - Point2f((0,radius))
        CairoMakie.arc!(ax,centre, radius, -1.2*pi,0.2*pi; arrow_attr...)
        CairoMakie.scatter!(ax,centre + radius*Point2f((cos(0.2*pi),sin(0.2*pi))), marker = arrow_type, color = arrow_attr.color, rotations = 0.2*pi, markersize = arrow_size)

        if ! isnothing(annotation)
            mid_point_text = centre - 1.5*Point2f((0,radius))
            CairoMakie.text!(ax,mid_point_text,text = annotation, align = (:center,:baseline),color = arrow_attr.color)
        end
    
    elseif node_id == :C
        centre = pos_start + Point2f((radius,0))
        CairoMakie.arc!(ax,centre, radius, -0.65*pi,0.65*pi; arrow_attr...)
        CairoMakie.scatter!(ax,centre + radius*Point2f((cos(0.65*pi),sin(0.65*pi))), marker = arrow_type, color = arrow_attr.color, rotations = 0.65*pi, markersize = arrow_size)

        if ! isnothing(annotation)
            mid_point_text = centre + 1.5*Point2f((0,radius))
            CairoMakie.text!(ax,mid_point_text,text = annotation, align = (:center,:baseline),color = arrow_attr.color)
        end

    elseif node_id == :A
        centre = pos_start - Point2f((radius,0))
        CairoMakie.arc!(ax,centre, radius, 0.35*pi,1.65*pi; arrow_attr...)
        CairoMakie.scatter!(ax,centre + radius*Point2f((cos(1.65*pi),sin(1.65*pi))), marker = arrow_type, color = arrow_attr.color, rotations = 1.65*pi, markersize = arrow_size)

        if ! isnothing(annotation)
            mid_point_text = centre - 1.5*Point2f((radius,0))
            CairoMakie.text!(ax,mid_point_text,text = annotation, align = (:center,:baseline),color = arrow_attr.color)
        end
    end

end

##########################


function draw_grn!(ax,network,draw_config::DrawGRNConfig,node_colors,fs,annotate = false,gene_letters = false)

    network_m = reshape(network,(3,4))

    non_zero_weights = Tuple.(findall(network_m .!= 0))

    if gene_letters
        CairoMakie.scatter!(ax,draw_config.fixed_layout_scatter,markersize = draw_config.node_size, color = :white, strokewidth = [draw_config.l,draw_config.l,draw_config.l,draw_config.l] , strokecolor = node_colors)
        CairoMakie.scatter!(ax,draw_config.fixed_layout_scatter,marker = ['A','B','C','M'], markersize = fs, color = node_colors)
    else
        CairoMakie.scatter!(ax,draw_config.fixed_layout_scatter,markersize = draw_config.node_size, color = node_colors)
    end

    symbolic_name = [:A,:B,:C,:M]

    annotations = string.(round.(network_m,digits = 2))

    for (pos_end_id,pos_start_id) in non_zero_weights

        reverse = (pos_start_id,pos_end_id) ∈ non_zero_weights

        pos_start = draw_config.fixed_layout_scatter[pos_start_id]
        pos_end = draw_config.fixed_layout_scatter[pos_end_id]

        pos_start_name = symbolic_name[pos_start_id]
        pos_end_name = symbolic_name[pos_end_id]

        if annotate
            if pos_start_name == pos_end_name
                if network_m[pos_end_id,pos_start_id] < 0
                    plot_self_arc!(ax,'_',draw_config.arrow_size_inh,draw_config.self_arc_radius,pos_start,annotations[pos_end_id,pos_start_id],pos_start_name,draw_config.arrow_attr)
                else
                    plot_self_arc!(ax,'▲',draw_config.arrow_size_act,draw_config.self_arc_radius,pos_start,annotations[pos_end_id,pos_start_id],pos_start_name,draw_config.arrow_attr)
                end
            else
                if pos_start_name == :M
                    if network_m[pos_end_id,pos_start_id] < 0
                        plot_shifted_arrows!(ax,'_',draw_config.arrow_size_inh,pos_start,pos_end,annotations[pos_end_id,pos_start_id],0.9*draw_config.shift_factor_node,draw_config.shift_factor_edge(false),-draw_config.shift_factor_label,draw_config.arrow_attr)
                    else
                        plot_shifted_arrows!(ax,'▲',draw_config.arrow_size_act,pos_start,pos_end,annotations[pos_end_id,pos_start_id],0.9*draw_config.shift_factor_node,draw_config.shift_factor_edge(false),-draw_config.shift_factor_label,draw_config.arrow_attr)
                    end
                else
                    if reverse
                        if network_m[pos_end_id,pos_start_id] < 0
                            plot_shifted_arrows!(ax,'_',draw_config.arrow_size_inh,pos_start,pos_end,annotations[pos_end_id,pos_start_id],draw_config.shift_factor_node,draw_config.shift_factor_edge(true),draw_config.shift_factor_label,draw_config.arrow_attr)
                        else
                            plot_shifted_arrows!(ax,'▲',draw_config.arrow_size_act,pos_start,pos_end,annotations[pos_end_id,pos_start_id],draw_config.shift_factor_node,draw_config.shift_factor_edge(true),draw_config.shift_factor_label,draw_config.arrow_attr)
                        end
                    else
                        if network_m[pos_end_id,pos_start_id] < 0
                            plot_shifted_arrows!(ax,'_',draw_config.arrow_size_inh,pos_start,pos_end,annotations[pos_end_id,pos_start_id],draw_config.shift_factor_node,draw_config.shift_factor_edge(false),draw_config.shift_factor_label,draw_config.arrow_attr)
                        else
                            plot_shifted_arrows!(ax,'▲',draw_config.arrow_size_act,pos_start,pos_end,annotations[pos_end_id,pos_start_id],draw_config.shift_factor_node,draw_config.shift_factor_edge(false),draw_config.shift_factor_label,draw_config.arrow_attr)
                        end
                    end
                end
            end
        else
            if pos_start_name == pos_end_name
                if network_m[pos_end_id,pos_start_id] < 0
                    plot_self_arc!(ax,'_',draw_config.arrow_size_inh,draw_config.self_arc_radius,pos_start,nothing,pos_start_name,draw_config.arrow_attr)
                else
                    plot_self_arc!(ax,'▲',draw_config.arrow_size_act,draw_config.self_arc_radius,pos_start,nothing,pos_start_name,draw_config.arrow_attr)
                end
            else
                if pos_start_name == :M
                    if network_m[pos_end_id,pos_start_id] < 0
                        plot_shifted_arrows!(ax,'_',draw_config.arrow_size_inh,pos_start,pos_end,nothing,0.9*draw_config.shift_factor_node,draw_config.shift_factor_edge(false),-draw_config.shift_factor_label,draw_config.arrow_attr)
                    else
                        plot_shifted_arrows!(ax,'▲',draw_config.arrow_size_act,pos_start,pos_end,nothing,0.9*draw_config.shift_factor_node,draw_config.shift_factor_edge(false),-draw_config.shift_factor_label,draw_config.arrow_attr)
                    end
                else
                    if reverse
                        if network_m[pos_end_id,pos_start_id] < 0
                            plot_shifted_arrows!(ax,'_',draw_config.arrow_size_inh,pos_start,pos_end,nothing,draw_config.shift_factor_node,draw_config.shift_factor_edge(true),draw_config.shift_factor_label,draw_config.arrow_attr)
                        else
                            plot_shifted_arrows!(ax,'▲',draw_config.arrow_size_act,pos_start,pos_end,nothing,draw_config.shift_factor_node,draw_config.shift_factor_edge(true),draw_config.shift_factor_label,draw_config.arrow_attr)
                        end
                    else
                        if network_m[pos_end_id,pos_start_id] < 0
                            plot_shifted_arrows!(ax,'_',draw_config.arrow_size_inh,pos_start,pos_end,nothing,draw_config.shift_factor_node,draw_config.shift_factor_edge(false),draw_config.shift_factor_label,draw_config.arrow_attr)
                        else
                            plot_shifted_arrows!(ax,'▲',draw_config.arrow_size_act,pos_start,pos_end,nothing,draw_config.shift_factor_node,draw_config.shift_factor_edge(false),draw_config.shift_factor_label,draw_config.arrow_attr)
                        end
                    end
                end
            end
        end

    end

    CairoMakie.xlims!(ax,draw_config.A_pos[1]-draw_config.shift_axis_edge,draw_config.C_pos[1]+draw_config.shift_axis_edge)
    CairoMakie.ylims!(ax,draw_config.B_pos[2]-draw_config.shift_axis_edge,draw_config.M_pos[2]+draw_config.shift_axis_edge_m)

    hidedecorations!(ax); hidespines!(ax)

end

function draw_grn!(ax,network,draw_config::DrawGRNConfig,node_colors,fs,annotations,annotate = true,gene_letters = false)

    network_m = reshape(network,(3,4))

    non_zero_weights = Tuple.(findall(network_m .!= 0))

    if gene_letters
        CairoMakie.scatter!(ax,draw_config.fixed_layout_scatter,markersize = draw_config.node_size, color = :white, strokewidth = [draw_config.l,draw_config.l,draw_config.l,draw_config.l] , strokecolor = node_colors)
        CairoMakie.scatter!(ax,draw_config.fixed_layout_scatter,marker = ['A','B','C','M'], markersize = fs, color = node_colors)
    else
        CairoMakie.scatter!(ax,draw_config.fixed_layout_scatter,markersize = draw_config.node_size, color = node_colors)
    end

    symbolic_name = [:A,:B,:C,:M]

    for (pos_end_id,pos_start_id) in non_zero_weights

        reverse = (pos_start_id,pos_end_id) ∈ non_zero_weights

        pos_start = draw_config.fixed_layout_scatter[pos_start_id]
        pos_end = draw_config.fixed_layout_scatter[pos_end_id]

        pos_start_name = symbolic_name[pos_start_id]
        pos_end_name = symbolic_name[pos_end_id]

        if annotate
            if pos_start_name == pos_end_name
                if network_m[pos_end_id,pos_start_id] < 0
                    plot_self_arc!(ax,'_',draw_config.arrow_size_inh,draw_config.self_arc_radius,pos_start,annotations[pos_end_id,pos_start_id],pos_start_name,draw_config.arrow_attr)
                else
                    plot_self_arc!(ax,'▲',draw_config.arrow_size_act,draw_config.self_arc_radius,pos_start,annotations[pos_end_id,pos_start_id],pos_start_name,draw_config.arrow_attr)
                end
            else
                if pos_start_name == :M
                    if network_m[pos_end_id,pos_start_id] < 0
                        plot_shifted_arrows!(ax,'_',draw_config.arrow_size_inh,pos_start,pos_end,annotations[pos_end_id,pos_start_id],0.9*draw_config.shift_factor_node,draw_config.shift_factor_edge(false),-draw_config.shift_factor_label,draw_config.arrow_attr)
                    else
                        plot_shifted_arrows!(ax,'▲',draw_config.arrow_size_act,pos_start,pos_end,annotations[pos_end_id,pos_start_id],0.9*draw_config.shift_factor_node,draw_config.shift_factor_edge(false),-draw_config.shift_factor_label,draw_config.arrow_attr)
                    end
                else
                    if reverse
                        if network_m[pos_end_id,pos_start_id] < 0
                            plot_shifted_arrows!(ax,'_',draw_config.arrow_size_inh,pos_start,pos_end,annotations[pos_end_id,pos_start_id],draw_config.shift_factor_node,draw_config.shift_factor_edge(true),draw_config.shift_factor_label,draw_config.arrow_attr)
                        else
                            plot_shifted_arrows!(ax,'▲',draw_config.arrow_size_act,pos_start,pos_end,annotations[pos_end_id,pos_start_id],draw_config.shift_factor_node,draw_config.shift_factor_edge(true),draw_config.shift_factor_label,draw_config.arrow_attr)
                        end
                    else
                        if network_m[pos_end_id,pos_start_id] < 0
                            plot_shifted_arrows!(ax,'_',draw_config.arrow_size_inh,pos_start,pos_end,annotations[pos_end_id,pos_start_id],draw_config.shift_factor_node,draw_config.shift_factor_edge(false),draw_config.shift_factor_label,draw_config.arrow_attr)
                        else
                            plot_shifted_arrows!(ax,'▲',draw_config.arrow_size_act,pos_start,pos_end,annotations[pos_end_id,pos_start_id],draw_config.shift_factor_node,draw_config.shift_factor_edge(false),draw_config.shift_factor_label,draw_config.arrow_attr)
                        end
                    end
                end
            end
        else
            if pos_start_name == pos_end_name
                if network_m[pos_end_id,pos_start_id] < 0
                    plot_self_arc!(ax,'_',draw_config.arrow_size_inh,draw_config.self_arc_radius,pos_start,nothing,pos_start_name,draw_config.arrow_attr)
                else
                    plot_self_arc!(ax,'▲',draw_config.arrow_size_act,draw_config.self_arc_radius,pos_start,nothing,pos_start_name,draw_config.arrow_attr)
                end
            else
                if pos_start_name == :M
                    if network_m[pos_end_id,pos_start_id] < 0
                        plot_shifted_arrows!(ax,'_',draw_config.arrow_size_inh,pos_start,pos_end,nothing,0.9*draw_config.shift_factor_node,draw_config.shift_factor_edge(false),-draw_config.shift_factor_label,draw_config.arrow_attr)
                    else
                        plot_shifted_arrows!(ax,'▲',draw_config.arrow_size_act,pos_start,pos_end,nothing,0.9*draw_config.shift_factor_node,draw_config.shift_factor_edge(false),-draw_config.shift_factor_label,draw_config.arrow_attr)
                    end
                else
                    if reverse
                        if network_m[pos_end_id,pos_start_id] < 0
                            plot_shifted_arrows!(ax,'_',draw_config.arrow_size_inh,pos_start,pos_end,nothing,draw_config.shift_factor_node,draw_config.shift_factor_edge(true),draw_config.shift_factor_label,draw_config.arrow_attr)
                        else
                            plot_shifted_arrows!(ax,'▲',draw_config.arrow_size_act,pos_start,pos_end,nothing,draw_config.shift_factor_node,draw_config.shift_factor_edge(true),draw_config.shift_factor_label,draw_config.arrow_attr)
                        end
                    else
                        if network_m[pos_end_id,pos_start_id] < 0
                            plot_shifted_arrows!(ax,'_',draw_config.arrow_size_inh,pos_start,pos_end,nothing,draw_config.shift_factor_node,draw_config.shift_factor_edge(false),draw_config.shift_factor_label,draw_config.arrow_attr)
                        else
                            plot_shifted_arrows!(ax,'▲',draw_config.arrow_size_act,pos_start,pos_end,nothing,draw_config.shift_factor_node,draw_config.shift_factor_edge(false),draw_config.shift_factor_label,draw_config.arrow_attr)
                        end
                    end
                end
            end
        end

    end

    CairoMakie.xlims!(ax,draw_config.A_pos[1]-draw_config.shift_axis_edge,draw_config.C_pos[1]+draw_config.shift_axis_edge)
    CairoMakie.ylims!(ax,draw_config.B_pos[2]-draw_config.shift_axis_edge,draw_config.M_pos[2]+draw_config.shift_axis_edge_m)

    hidedecorations!(ax); hidespines!(ax)

end


function draw_grn!(ax,network,reverse_list,draw_config::DrawGRNConfig,node_colors,fs,annotate = false,gene_letters = false)

    network_m = reshape(network,(3,4))

    non_zero_weights = Tuple.(findall(network_m .!= 0))

    if gene_letters
        CairoMakie.scatter!(ax,draw_config.fixed_layout_scatter,markersize = draw_config.node_size, color = :white, strokewidth = [draw_config.l,draw_config.l,draw_config.l,draw_config.l] , strokecolor = node_colors)
        CairoMakie.scatter!(ax,draw_config.fixed_layout_scatter,marker = ['A','B','C','M'], markersize = fs, color = node_colors)
    else
        CairoMakie.scatter!(ax,draw_config.fixed_layout_scatter,markersize = draw_config.node_size, color = node_colors)
    end

    symbolic_name = [:A,:B,:C,:M]

    annotations = string.(round.(network_m,digits = 2))

    for (pos_end_id,pos_start_id) in non_zero_weights

        reverse = (pos_start_id,pos_end_id) ∈ reverse_list

        pos_start = draw_config.fixed_layout_scatter[pos_start_id]
        pos_end = draw_config.fixed_layout_scatter[pos_end_id]

        pos_start_name = symbolic_name[pos_start_id]
        pos_end_name = symbolic_name[pos_end_id]

        if annotate
            if pos_start_name == pos_end_name
                if network_m[pos_end_id,pos_start_id] < 0
                    plot_self_arc!(ax,'_',draw_config.arrow_size_inh,draw_config.self_arc_radius,pos_start,annotations[pos_end_id,pos_start_id],pos_start_name,draw_config.arrow_attr)
                else
                    plot_self_arc!(ax,'▲',draw_config.arrow_size_act,draw_config.self_arc_radius,pos_start,annotations[pos_end_id,pos_start_id],pos_start_name,draw_config.arrow_attr)
                end
            else
                if pos_start_name == :M
                    if network_m[pos_end_id,pos_start_id] < 0
                        plot_shifted_arrows!(ax,'_',draw_config.arrow_size_inh,pos_start,pos_end,annotations[pos_end_id,pos_start_id],0.9*draw_config.shift_factor_node,draw_config.shift_factor_edge(false),-draw_config.shift_factor_label,draw_config.arrow_attr)
                    else
                        plot_shifted_arrows!(ax,'▲',draw_config.arrow_size_act,pos_start,pos_end,annotations[pos_end_id,pos_start_id],0.9*draw_config.shift_factor_node,draw_config.shift_factor_edge(false),-draw_config.shift_factor_label,draw_config.arrow_attr)
                    end
                else
                    if reverse
                        if network_m[pos_end_id,pos_start_id] < 0
                            plot_shifted_arrows!(ax,'_',draw_config.arrow_size_inh,pos_start,pos_end,annotations[pos_end_id,pos_start_id],draw_config.shift_factor_node,draw_config.shift_factor_edge(true),draw_config.shift_factor_label,draw_config.arrow_attr)
                        else
                            plot_shifted_arrows!(ax,'▲',draw_config.arrow_size_act,pos_start,pos_end,annotations[pos_end_id,pos_start_id],draw_config.shift_factor_node,draw_config.shift_factor_edge(true),draw_config.shift_factor_label,draw_config.arrow_attr)
                        end
                    else
                        if network_m[pos_end_id,pos_start_id] < 0
                            plot_shifted_arrows!(ax,'_',draw_config.arrow_size_inh,pos_start,pos_end,annotations[pos_end_id,pos_start_id],draw_config.shift_factor_node,draw_config.shift_factor_edge(false),draw_config.shift_factor_label,draw_config.arrow_attr)
                        else
                            plot_shifted_arrows!(ax,'▲',draw_config.arrow_size_act,pos_start,pos_end,annotations[pos_end_id,pos_start_id],draw_config.shift_factor_node,draw_config.shift_factor_edge(false),draw_config.shift_factor_label,draw_config.arrow_attr)
                        end
                    end
                end
            end
        else
            if pos_start_name == pos_end_name
                if network_m[pos_end_id,pos_start_id] < 0
                    plot_self_arc!(ax,'_',draw_config.arrow_size_inh,draw_config.self_arc_radius,pos_start,nothing,pos_start_name,draw_config.arrow_attr)
                else
                    plot_self_arc!(ax,'▲',draw_config.arrow_size_act,draw_config.self_arc_radius,pos_start,nothing,pos_start_name,draw_config.arrow_attr)
                end
            else
                if pos_start_name == :M
                    if network_m[pos_end_id,pos_start_id] < 0
                        plot_shifted_arrows!(ax,'_',draw_config.arrow_size_inh,pos_start,pos_end,nothing,0.9*draw_config.shift_factor_node,draw_config.shift_factor_edge(false),-draw_config.shift_factor_label,draw_config.arrow_attr)
                    else
                        plot_shifted_arrows!(ax,'▲',draw_config.arrow_size_act,pos_start,pos_end,nothing,0.9*draw_config.shift_factor_node,draw_config.shift_factor_edge(false),-draw_config.shift_factor_label,draw_config.arrow_attr)
                    end
                else
                    if reverse
                        if network_m[pos_end_id,pos_start_id] < 0
                            plot_shifted_arrows!(ax,'_',draw_config.arrow_size_inh,pos_start,pos_end,nothing,draw_config.shift_factor_node,draw_config.shift_factor_edge(true),draw_config.shift_factor_label,draw_config.arrow_attr)
                        else
                            plot_shifted_arrows!(ax,'▲',draw_config.arrow_size_act,pos_start,pos_end,nothing,draw_config.shift_factor_node,draw_config.shift_factor_edge(true),draw_config.shift_factor_label,draw_config.arrow_attr)
                        end
                    else
                        if network_m[pos_end_id,pos_start_id] < 0
                            plot_shifted_arrows!(ax,'_',draw_config.arrow_size_inh,pos_start,pos_end,nothing,draw_config.shift_factor_node,draw_config.shift_factor_edge(false),draw_config.shift_factor_label,draw_config.arrow_attr)
                        else
                            plot_shifted_arrows!(ax,'▲',draw_config.arrow_size_act,pos_start,pos_end,nothing,draw_config.shift_factor_node,draw_config.shift_factor_edge(false),draw_config.shift_factor_label,draw_config.arrow_attr)
                        end
                    end
                end
            end
        end

    end

    CairoMakie.xlims!(ax,draw_config.A_pos[1]-draw_config.shift_axis_edge,draw_config.C_pos[1]+draw_config.shift_axis_edge)
    CairoMakie.ylims!(ax,draw_config.B_pos[2]-draw_config.shift_axis_edge,draw_config.M_pos[2]+draw_config.shift_axis_edge_m)

    hidedecorations!(ax); hidespines!(ax)

end


function draw_grn_mutant!(ax,network_old,network_new,draw_config::DrawGRNConfig,draw_config_generator,node_colors,fs,annotate = false,gene_letters = false)

    unchanged_network = Int.(network_old .== network_new) .* network_new
    network_additions = Int.(network_old .!= network_new) .* network_new

    network_new_m = reshape(network_new,(3,4))

    reverse_list = Tuple.(findall(network_new_m .!= 0))

    draw_grn!(ax,unchanged_network,reverse_list,draw_config,node_colors,fs,annotate,gene_letters)

    arrow_attr_new = (; linewidth = draw_config.arrow_attr.linewidth, color = :red)

    draw_config_new = draw_config_generator()

    draw_config_new.arrow_attr = arrow_attr_new

    draw_grn!(ax,network_additions,reverse_list,draw_config_new,node_colors,fs,annotate,gene_letters)

end