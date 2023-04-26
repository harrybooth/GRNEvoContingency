# MOL: u_{j}(t) = u(x_j,t) where x_j = j*dx

function gene_regulation_1d!(dg,g,p,t)

    w, λg = p

    @inbounds for j in 1:Nc
        x = tissue[j]
        dg[1,j] = σ(w[1,1]*g[1,j]+ w[1,2]*g[2,j] + w[1,3]*g[3,j] + w[1,4]*morph(x)) - λg[1]*g[1,j]
        dg[2,j] = σ(w[2,1]*g[1,j]+ w[2,2]*g[2,j] + w[2,3]*g[3,j] + w[2,4]*morph(x)) - λg[2]*g[2,j]
        dg[3,j] = σ(w[3,1]*g[1,j]+ w[3,2]*g[2,j] + w[3,3]*g[3,j] + w[3,4]*morph(x)) - λg[3]*g[3,j]
    end
    
end

function init_gene_regulation_1d(start_conc)
    return start_conc .* ones(Ng,Nc)
end

# For homotopy methods

s(i,w,g,m) = sum(w[i,k] * g[k] for k in 1:3) + w[i, 4] * m

δ(i,w,g,m,λg) = σ(s(i,w,g,m)) - λg*g[i]

function Δ(i,w,g,d,m,λg)
    I = s(i,w,g,m) + h_a
    [
      # d[j,i] = sqrt(I^2 + 1)
      d[i]^2 - (I^2 + h_b)
      # σ(s(i,j)) - λ_g g_j^{(i)}
      0.5 * ((I / d[i]) + 1) - λg * g[i]
    ]
  end
  
function GRN(w,g,d,m,λg)
    System(
        [Δ(1,w,g,d,m,λg)
         Δ(2,w,g,d,m,λg)
         Δ(3,w,g,d,m,λg)];
      variables=[vec(g);vec(d)],
      parameters=[vec(w);λg;m])
end

function GRN_a(w,g,d,λg)
    System(reduce(vcat,[
        [Δ(1,w,g,d,morph(tissue[cell]),λg)
         Δ(2,w,g,d,morph(tissue[cell]),λg)
         Δ(3,w,g,d,morph(tissue[cell]),λg)] for cell in 1:Nc
      ]);
      variables=[vec(g);vec(d)],
      parameters=[vec(w);λg])
  end
  




