using HomotopyContinuation

m = 4
N = 1


@var w[1:3, 1:m] M[1:N] λ_g g[1:3, 1:N]

@var d[1:3, 1:N]


s(i, j) = sum(w[j, k] * g[k, i] for k in 1:3) + w[j, m] * M[i]


function Δ(i, j)
  I = s(i, j)
  [
    # d[j,i] = sqrt(I^2 + 1)
    d[j, i]^2 - (I^2 + 1)
    # σ(s(i,j)) - λ_g g_j^{(i)}
    0.5 * ((I / d[j, i]) + 1) - λ_g * g[j, i]
  ]
end

function system(i)
  System([
      Δ(i, 1)
      Δ(i, 2)
      Δ(i, 3)
    ];
    variables=[g[:, i]; d[:, i]],
    parameters=[vec(w); vec(M); λ_g])
end

F = system(1)

r = monodromy_solve(F)
S = solutions(r)
q = parameters(r)


# Choose something sensible here
P = [randn(14) for k in 1:10]

# Track to all the different parameters
P_results = solve(F, S; start_parameters=q, target_parameters=P)

# Comment / Hint
# monodromy_solve choose random complex values for all parameters (from a normal distribution) (this is `q` above)
# If your real parameter values have a different scale, then you should probably move first the generic solutions `S`
# to a new set of random complex parameters of the same scale as your real parameters.
