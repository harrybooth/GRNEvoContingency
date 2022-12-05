using HomotopyContinuation

const Nc = 50

const tissue = range(0,1,length = Nc)

const dx = step(tissue)

const  λ_g = 0.05

const h_a = -1.
const h_b = 0.1

σ(I) = 0.5*(((I + h_a)/sqrt((I + h_a)^2+h_b)) + 1)

morph(x) = 10*exp(-x/0.4)

@var w[1:3, 1:4]  g[1:3,1:Nc], d[1:3,1:Nc]

s(i, cell) = sum(w[i,k] * g[k, cell] for k in 1:3) + w[i, 4] * morph(tissue[cell])

function Δ(i,cell)
  I = s(i,cell) + h_a
  [
    # d[j,i] = sqrt(I^2 + 1)
    d[i,cell]^2 - (I^2 + h_b)
    # σ(s(i,j)) - λ_g g_j^{(i)}
    0.5 * ((I / d[i,cell]) + 1) - λ_g * g[i,cell]
  ]
end

function GRN()
  System(reduce(vcat,[
      [Δ(1,cell)
      Δ(2,cell)
      Δ(3,cell)] for cell in 1:Nc
    ]);
    variables=[vec(g); vec(d)],
    parameters=[vec(w);])
end

F = GRN()

w_ex = randn(ComplexF64, length(vec(w)));

result_start = solve(F, target_parameters = w_ex)