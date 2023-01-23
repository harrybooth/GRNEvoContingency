using PyCall

fred = pyimport("Fred")

################### Utility functions

function refresh_type(v)
    [i for i in v]
end

function attach_column_to_tuples(column::Vector{Float64},tuples::Vector{Tuple{Float64,Float64}})
    return hcat(column,hcat(first.(tuples), last.(tuples)))
end

function tuples_to_matrix(tuples::Vector{Tuple{Float64,Float64}})
    return hcat(first.(tuples), last.(tuples))
end

################### Frechet distances

function frechet(metric::PreMetric, P::AbstractMatrix, Q::AbstractMatrix)
    dP, m = size(P)
    dQ, n = size(Q)

    dP != dQ && throw(DimensionMismatch(
        "Points in polygonal lines P and Q must have the same dimension."))

    couplings = pairwise(metric, P, Q, dims=2)

    @inbounds for i in 2:m
        couplings[i, 1] = max(couplings[i-1, 1], couplings[i, 1])
    end

    @inbounds for j in 2:n
        couplings[1, j] = max(couplings[1, j-1], couplings[1, j])
        for i in 2:m
            carried_coupling = min(couplings[i-1, j-1],
                                   couplings[i, j-1],
                                   couplings[i-1, j])
            couplings[i, j] = max(carried_coupling, couplings[i, j])
        end
    end

    return couplings[m, n]
end

function pairwise_frechet(metric::PreMetric,X)

    result = zeros(length(X),length(X))

    for j in 1:size(result,2)
        for i in 1:size(result,1)
            if i < j
                result[i,j] = 1.
            end
        end
    end

    choices = Tuple.(findall(x->x>0,result))

    for entry in choices
        result[entry...] = frechet(metric,X[entry[1]],X[entry[2]])
    end
    
    return Symmetric(result) |> collect
end

function pairwise_frechet(metric::PreMetric,X,Y)

    result = zeros(length(X),length(Y))

    for j in 1:size(result,2)
        for i in 1:size(result,1)
            if i <= j
                result[i,j] = 1.
            end
        end
    end

    choices = Tuple.(findall(x->x>0,result))

    for entry in choices
        result[entry...] = frechet(metric,X[entry[1]],Y[entry[2]])
    end
    
    return Symmetric(result) |> collect
end

function pairwise_cts_frechet(X)

    result = zeros(length(X),length(X))

    X_t = map(x->x |> transpose |> collect, X);

    for j in 1:size(result,2)
        for i in 1:size(result,1)
            if i <= j
                result[i,j] = 1.
            end
        end
    end

    choices = Tuple.(findall(x->x>0,result))

    for entry in choices
        result[entry...] = fred.continuous_frechet(fred.Curve(X_t[entry[1]]), fred.Curve(X_t[entry[2]])).value
    end
    
    return Symmetric(result) |> collect
end

################### Distance functions

function compute_distance_mat!(Tr::Union{PhenoTrajectories,GenoTrajectories},metric::PreMetric)
    Tr.distance_mat  = pairwise_frechet(metric,Tr.traj_m)
end

function compute_cts_distance_mat!(Tr::Union{PhenoTrajectories,GenoTrajectories},metric::PreMetric)
    Tr.distance_mat  = pairwise_cts_frechet(Tr.traj_m)
end

################### Performance measures

function calinski_harabasz_dmat(distance_mat,assignments)

    n_clust = unique(assignments) |> length
    n_sample = length(assignments)

    trajectory_exemplar = argmin(sum(distance_mat, dims = 2))[1]

    cluster_exemplars = zeros(Int64,n_clust)

    bgss = 0.
    wgss = 0.

    cluster_exemplars = zeros(n_clust)

    for cl in unique(assignments)

        n_cl = sum(assignments .== cl)

        pt_cl = copy(distance_mat)
        pt_cl[:,assignments .!= cl] .= 0 # no contribution to row sum from points outside cluster
        pt_cl[assignments .!= cl,:] .= Inf # rows outside of the cluster can never be selected as centroid
        cl_exemplar = argmin(sum(pt_cl, dims = 2))[1] # From modified distance matrix, find row with the smallest summed distances, hence corresponding to the most "central" trajectory

        cluster_exemplars[cl] = cl_exemplar
        
        bgss += n_cl*distance_mat[trajectory_exemplar,cl_exemplar]

        wgss += sum(distance_mat[assignments .==  cl,cl_exemplar]) 

    end

    return (bgss / wgss) * ((n_sample - n_clust)/(n_clust - 1))

end


function between_like_within_like(distance_mat,assignments)

    n_clust = unique(assignments) |> length
    n_sample = length(assignments)

    trajectory_exemplar = argmin(sum(distance_mat, dims = 2))[1]

    cluster_exemplars = zeros(Int64,n_clust)

    WL = 0
    BL = 0

    for cl in unique(assignments)

        n_cl = sum(assignments .== cl)

        pt_cl = copy(distance_mat)
        pt_cl[:,assignments .!= cl] .= 0 # no contribution to row sum from points outside cluster
        pt_cl[assignments .!= cl,:] .= Inf # rows outside of the cluster can never be selected as centroid
        cl_exemplar = argmin(sum(pt_cl, dims = 2))[1] # From modified distance matrix, find row with the smallest summed distances, hence corresponding to the most "central" trajectory

        cluster_exemplars[cl] = cl_exemplar

        WL += sum(distance_mat[assignments .==  cl,cl_exemplar]) / n_cl
        BL += distance_mat[trajectory_exemplar,cl_exemplar]

    end


    return WL / n_clust, BL / n_clust , cluster_exemplars, trajectory_exemplar

end

################### Plotting

function return_arrows(traj::Vector{Tuple{Float64,Float64}})

    traj_m = tuples_to_matrix(traj)

    all_diff_x = []
    all_diff_y = []

    x_s = []
    y_s = []

    for di in 1:size(traj_m,1)-1
        diff_x = traj_m[di+1,1] - traj_m[di,1]
        diff_y = traj_m[di+1,2] - traj_m[di,2]
        push!(all_diff_x,diff_x)
        push!(all_diff_y,diff_y)
        push!(x_s,traj_m[di,1])
        push!(y_s,traj_m[di,2])

    end

    return x_s, y_s, all_diff_x,all_diff_y
end

function cluster_density(pt :: Union{PhenoTrajectories,GenoTrajectories}, n_clusters::Int64)

    cm = countmap(pt.cluster_assignments)
    result = zeros(n_clusters)
    for cl in 1:n_clusters
        result[cl] = cm[cl] / length(pt.cluster_assignments)
    end
    return result
end

function plot_methodology_pheno_traj(pt::PhenoTrajectories, sample_choices::Vector{Int64}, topology_name::String, target = (40,20))

    fig = CairoMakie.Figure(resolution = (800, 800),fontsize = 27.)

    # sample_choices = [12,8]
    # sample_choices = [2,82]
    
    color_choices = palette(:tab10)

    # Contour plot

    xmin, x_max = (0,120)
    ymin, y_max = (0,120)

    c = LinRange(xmin, x_max, 200)
    w = LinRange(ymin, y_max, 200)

    z = -1 .* ( (c' .- target[1]).^2 .+ (w .- target[2]).^2 )

    gl = fig[1,1]

    ax1  = Axis(gl[1,1], backgroundcolor = "white", xlabel = "Centre position (cell id)", ylabel = "Stripe width (cells)", title = "Example phenotype trajectories for " * topology_name * " topology")

    c = CairoMakie.contourf!(ax1,c,w,transpose(z),levels = 10)

    path_n = 1

    for (sample_id,cc) in zip(sample_choices,color_choices)
        CairoMakie.scatter!(ax1,unique(pt.traj[sample_id]), label = "Trajectory " * string(path_n),markersize = 15., marker = :x,color = cc)
        CairoMakie.arrows!(ax1,return_arrows(unique(pt.traj[sample_id]))...,arrowsize  = 8.)

        path_n +=1
    end

    CairoMakie.scatter!(ax1,target,markersize = 20., marker = :●,label = "Target Phenotype", color = :purple)

    CairoMakie.axislegend()

    Colorbar(gl[1,2], c)

    fig
end
