using HomotopyContinuation
using Symbolics
using LinearAlgebra
using Plots

const Nc = 200

const tissue = range(0,1,length = Nc)

const dx = step(tissue)

const  λ_g = 0.05

const h_a = -1.
const h_b = 0.1

σ(I) = 0.5*(((I + h_a)/sqrt((I + h_a)^2+h_b)) + 1)

morph(x) = 10*exp(-x/0.4)

# Defining system for Homotopy continuation

@var w[1:3, 1:4]  g[1:3], d[1:3], m, λg

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
  

# Symbolics - create jacobian

λg_default = 0.05

@variables G[1:3],M,λG

@variables W[1:3,1:4]

J = Symbolics.jacobian([δ(1,W,G,M,λG),δ(2,W,G,M,λG),δ(3,W,G,M,λG)], G)

J_call = eval(Symbolics.build_function(J,[G;vec(W); λG; M])[1])

# Define functions which will specify whether a solution is valid or not (positive real +  all real eigenvalues negative)

function satisfies_system(G,p)

    G1 = G[1]
    G2 = G[2]
    G3 = G[3]

    w11, w21, w31, w12, w22, w32, w13, w23, w33, w1m, w2m, w3m, λg, M = p

    r1 = σ(w11*G1 + w12*G2 + w13*G3 + w1m*M) - λg*G1
    r2 = σ(w21*G1 + w22*G2 + w23*G3 + w2m*M) - λg*G2
    r3 = σ(w31*G1 + w32*G2 + w33*G3 + w3m*M) - λg*G3

    eig = J_call([G1,G2,G3,w11, w21, w31, w12, w22, w32, w13, w23, w33, w1m, w2m, w3m, λg,M])

    return all(abs.([r1,r2,r3]) .< 1e-11) && all(real.(eigvals(eig)) .< 0)
end

function get_valid_solutions(r,p)
    rp = filter(s -> all(s .> 0), real_solutions(r))
    rpv = filter(x->satisfies_system(x,p),rp)
    return rpv
end

# solve 1 cell

F = GRN(w,g,d,m,λg)

r = monodromy_solve(F)

S = solutions(r) # this seems to randomly determine the solutions I find in the methodology below
q = parameters(r);

# extend to Nc=200 cells using parameter homotopy, keeping real and stable solutions

w_ex =   [0.0 0.0 0.0 0.114353;
-0.560775 0.359711 0.0 0.0;
-0.0725802 -0.407225 0.527974 0.0]

w_ex = example_networks["overlap_dom"]

λg_ex = 0.05

weights = w_ex[:];

P = [[weights...,λg_ex,morph(x)] for x in tissue];

data_points = solve(F, S; start_parameters=q, target_parameters=P,transform_result = (r,p) ->get_valid_solutions(r,p));

d1 = map(y->y[1],filter(x->length(x) > 0,data_points))

plot(map(x->x[1],d1),label = "G1",title = "Alternative Method - classical")
plot!(map(x->x[2],d1),label = "G2")
plot!(map(x->x[3],d1),label = "G3")