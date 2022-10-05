function solve_SS_prob(prob_ode)
    solve(prob_ode,AutoTsit5(Rosenbrock23()), isoutofdomain=(u,p,t) -> any(x -> x < 0, u), callback = TerminateSteadyState(1e-5,1e-3),maxiters = 1e4,verbose = false, save_everystep = false);
end

function init_SS_prob(p)

    g0 = init_gene_regulation_1d(0.01)

    tspan = (0,Inf)

    dg0 = copy(g0) # what does dg look like

    jac_sparsity = Symbolics.jacobian_sparsity((dg,g)->gene_regulation_1d!(dg,g,p,0.0),dg0,g0)

    f = ODEFunction(gene_regulation_1d!;jac_prototype=float.(jac_sparsity))

    ODEProblem(f,g0,tspan,p)
end

function remake_connections(prob,w,p_old)
    remake(prob,p = (w,p_old[2:end]...))
end

pd = (u,p,t) -> any(x -> x < 0, u)
stab_condition = TerminateSteadyState(1e-5,1e-3);
