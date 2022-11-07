const Nc = 200
const Ng = 3
const L = 1.
const θ = 5.
const c0 = 10.
const λm = 0.4

# const h_a = 1
# const h_b = 0.1

const h_a = 1.
const h_b = 0.1

const tissue = range(0,L,length = Nc)

m(x) = c0*exp(-x/λm)

# σ(I) = 1/(1+exp(θ-θ*I))  # σ(0.) > 0 ?

# σ(I) = 0.5*(((I - h_a)/sqrt((I - h_a)^2+h_b)) + 1) # σ(0.) > 0 ?

σ(I) = 0.5*((I/sqrt(I^2+1)) + 1)

# MOL: u_{j}(t) = u(x_j,t) where x_j = j*dx

function gene_regulation_1d!(dg,g,p,t)

    w, λg = p

    @inbounds for j in 1:Nc
        x = tissue[j]
        dg[1,j] = σ(w[1,1]*g[1,j]+ w[1,2]*g[2,j] + w[1,3]*g[3,j] + w[1,4]*m(x)) - λg[1]*g[1,j]
        dg[2,j] = σ(w[2,1]*g[1,j]+ w[2,2]*g[2,j] + w[2,3]*g[3,j] + w[2,4]*m(x)) - λg[2]*g[2,j]
        dg[3,j] = σ(w[3,1]*g[1,j]+ w[3,2]*g[2,j] + w[3,3]*g[3,j] + w[3,4]*m(x)) - λg[3]*g[3,j]
    end
    
end


function init_gene_regulation_1d(start_conc)
    return start_conc .* ones(Ng,Nc)
end



