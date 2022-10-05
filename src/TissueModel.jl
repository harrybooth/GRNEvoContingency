using DifferentialEquations
# using Symbolics

const Nc = 100
const Ng = 3
const L = 1.
const θ = 5.
const c0 = 1.
const λm = 0.4

const tissue = range(0,L,length = Nc)

limit(a, N) = a == N+1 ? N : a == 0 ? 1 : a

m(x) = c0*exp(-x/λm)

σ(I) = 1/(1+exp(θ-θ*I))  # σ(0.) > 0 ?

# MOL: u_{j}(t) = u(x_j,t) where x_j = j*dx

function gene_regulation_1d!(dg,g,p,t)

    w, Dg, λg, dx = p

    h = 1/dx^2

    @inbounds for j in 1:Nc
        x = tissue[j]
        j_left, j_right = limit(j-1,Nc), limit(j+1,Nc)
        dg[1,j] = σ(w[1,1]*g[1,j]+ w[1,2]*g[2,j] + w[1,3]*g[3,j] + w[1,4]*m(x)) + h*Dg[1]*(g[1,j_left] + g[1,j_right] - 2*g[1,j]) - λg[1]*g[1,j]
        dg[2,j] = σ(w[2,1]*g[1,j]+ w[2,2]*g[2,j] + w[2,3]*g[3,j] + w[2,4]*m(x)) + h*Dg[2]*(g[2,j_left] + g[2,j_right] - 2*g[2,j]) - λg[2]*g[2,j]
        dg[3,j] = σ(w[3,1]*g[1,j]+ w[3,2]*g[2,j] + w[3,3]*g[3,j] + w[3,4]*m(x)) + h*Dg[3]*(g[3,j_left] + g[3,j_right] - 2*g[3,j]) - λg[3]*g[3,j]
    end
    
end

function init_gene_regulation_1d(start_conc)

    g0 = zeros(Ng,Nc)

    for j in 1:Nc
        for i in 1:Ng
            g0[i,j] = start_conc
        end
    end

    return g0
end


# const w_ex = [0.0 0.0 0.0 1.04715;
#     1.26868 0.0 0.0 0.0;
#     0.0370965 -0.281208 0.972401 0.0]

const w_ex = [ 0.0 0.0 0.0 1.25937;
0.167521 0.0 0.0 0.0;
0.257394 -2.84893 2.34484 0.0]

const Dg_ex = zeros(Ng)
const λg_ex = 0.05 .* ones(Ng)

