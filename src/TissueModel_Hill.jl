const Nc = 100
const Ng = 3
const L = 1.
const θ = 5.
const c0 = 10.
const λm = 0.4
const b = 1.

const tissue = range(0,L,length = Nc)

m(x) = c0*exp(-x/λm)

h2(w,g) = w >= 0 ? 1 + (w*g)^2/(b^2+(w*g)^2) : b^2/(b^2+(w*g)^2) # squared so dont need abs value of w

# MOL: u_{j}(t) = u(x_j,t) where x_j = j*dx

function gene_regulation_1d!(dg,g,p,t)

    w, λg = p

    @inbounds for j in 1:Nc
        x = tissue[j]
        dg[1,j] = h2(w[1,1],g[1,j])*h2(w[1,2],g[2,j])*h2(w[1,3],g[3,j])*h2(w[1,4],m(x)) - λg[1]*g[1,j]
        dg[2,j] = h2(w[2,1],g[1,j])*h2(w[2,2],g[2,j])*h2(w[2,3],g[3,j])*h2(w[2,4],m(x)) - λg[2]*g[2,j]
        dg[3,j] = h2(w[3,1],g[1,j])*h2(w[3,2],g[2,j])*h2(w[3,3],g[3,j])*h2(w[3,4],m(x)) - λg[3]*g[3,j]
    end
    
end

function init_gene_regulation_1d(start_conc)
    return start_conc .* ones(Ng,Nc)
end



