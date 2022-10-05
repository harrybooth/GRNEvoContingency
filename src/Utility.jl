using ColorTypes
using Printf
using BenchmarkTools

function visualize_steady_state(sol_ss_u)
    norm_ss = (sol_ss_u .- minimum(sol_ss_u)) ./ (maximum(sol_ss_u) .- minimum(sol_ss_u)) 
    mapslices(x->RGB{Float64}(x...),norm_ss,dims = 1)
end 

function plot_SS(sol_u_end)
    plot(transpose(sol_u_end),xlabel = "Tissue",ylabel = "Concentration",label = ["Gene A" "Gene B" "Gene C (output)"])
    plot!(map(x->m(x),tissue),label = "Morphogen (steady state)")
end