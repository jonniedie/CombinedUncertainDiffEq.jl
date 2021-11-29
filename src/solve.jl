"""
    uncomponentarrayize(prob; u0=prob.u0, p=prob.p)

Convert an `ODEProblem` whose function requires `ComponentArray`s for states and parameters
to one that only requires plain `Vector`s and converts them internally. This is needed
for `koopman_expectation`.
"""
function uncomponentarrayize(prob; u0=prob.u0, p=prob.p)
    f!  = prob.f.f
    ax_x = getaxes(u0)
    ax_p = getaxes(p)
    new_f! = function (dx, x, p, t; kwargs...)
        f!(ComponentArray(dx, ax_x), ComponentArray(x, ax_x), ComponentArray(p, ax_p), t; kwargs...)
    end
    # # TODO: Handle integrator, affect!, and affect_neg!
    # kwargs = if :callback in keys(prob.kwargs)
    #     cb = prob.kwargs[:callback]
    #     old_condition = deepcopy(cb.condition)
    #     condition = (x,t,integrator) -> old_condition(ComponentArray(x, ax_x), t, integrator)
    #     @set! cb.condition = condition
    #     (; prob.kwargs..., callback=cb)
    # else
    #     prob.kwargs
    # end
    kwargs = prob.kwargs
    return ODEProblem(new_f!, getdata(u0), prob.tspan, getdata(p); kwargs...)
end

componentarray_like(a, b::ComponentArray) = ComponentArray(a, getaxes(b))
componentarray_like(a::ComponentArray, b::ComponentArray) = a

"""
    mc_solve(prob, u0=prob.u0, p=prob.p, args...; trajectories, kwargs...)

Solve an `ODEProblem` `prob` with uncertain parameters and initial conditions as a Monte Carlo
simulation
"""
function mc_solve(prob, u0=prob.u0, p=prob.p, args...;
                    trajectories,
                    u0_CoV = (_u,_p) -> componentarray_like(_u, u0),
                    p_CoV = (_u,_p) -> componentarray_like(_p, p),
                    kwargs...)

    u0 = recursive_convert(to_distribution, u0)
    p = recursive_convert(to_distribution, p)
    prob_func = function (prob, i, repeat)
        _u0 = u0_CoV(_rand.(u0), p)
        _p = p_CoV(u0, _rand.(p))
        remake(prob, u0=_u0, p=_p)
    end
    mc_prob = EnsembleProblem(prob; prob_func=prob_func)
    return solve(mc_prob, args...; trajectories=trajectories, kwargs...)
end


"""
    mc_expectation(g, prob, u0, p, args...; trajectories, kwargs...)

Calculate the expectation of observable function `g` on the results of an `ODEProblem` `prob`
with uncertain initial conditions `u0` and `p`. This is the same as calling `expectation`
from DiffEqUncertainty.jl with the `Koopman` expectation algorithm with an added benefit of
being able to handle `ComponentArray` initial conditions and parameters as well as the
ability to specify either `Interval`s or `Distribution`s for uncertainty.
"""
function mc_expectation(g, prob, u0, p, args...; trajectories, kwargs...)
    sols = mc_solve(prob, u0, p, args...; trajectories, kwargs...)
    return mean(g, sols.u)
end


"""
    koopman_expectation(g, prob, u0, p, args...; kwargs...)

Calculate the expectation of observable function `g` on the results of an `ODEProblem` `prob`
with uncertain initial conditions `u0` and `p`. This is the same as calling `expectation`
from DiffEqUncertainty.jl with the `Koopman` expectation algorithm with an added benefit of
being able to handle `ComponentArray` initial conditions and parameters as well as the
ability to specify either `Interval`s or `Distribution`s for uncertainty.
"""
function koopman_expectation(g, prob, u0, p, args...; kwargs...)
    u0 = recursive_convert(to_distribution, u0)
    p = recursive_convert(to_distribution, p)
    return expectation(g, prob, u0, p, Koopman(), args...; kwargs...)
end

function koopman_expectation(g, prob, u0::ComponentArray, p::ComponentArray, args...; kwargs...)
    # new_prob = uncomponentarrayize(prob; u0=u0, p=p)
    ax_u0 = getaxes(u0)
    ax_p = getaxes(p)
    u0 = recursive_convert(to_distribution, getdata(u0))
    p = recursive_convert(to_distribution, getdata(p))
    u0_CoV = (u,p) -> ComponentArray(u, ax_u0)
    p_CoV = (u,p) -> ComponentArray(p, ax_p)
    return expectation(g, prob, u0, p, Koopman(), args...; u0_CoV, p_CoV, kwargs...)
end


"""
    reachability_solve(prob, u0, p, args...; kwargs...)

Solve an `ODEProblem` `prob` with uncertain parameters and initial conditions as a reachability
simulation. The output zonotope will enclose the entire reachable state.
"""
function reachability_solve(prob::ODEProblem, args...; kwargs...)
    new_prob = reachability_problem(prob)
    return RA.solve(new_prob, args...; prob.tspan, kwargs...)
end
reachability_solve(ivp::RA.InitialValueProblem, args...; kwargs...) = RA.solve(ivp, args...; kwargs...)
reachability_solve(f::Function, u0, tspan, p, args...; kwargs...) = reachability_solve(ODEProblem(f, u0, tspan, p), args...; kwargs...)


reachability_problem(prob::ODEProblem, args...; kwargs...) = reachability_problem(prob.f.f, prob.u0, prob.p, args...; inplace=isinplace(prob), kwargs...)
function reachability_problem(prob_f::Function, u0, p, args...; inplace=true, kwargs...)
    x0, unflatten_u, unflatten_p, unflatten_t = reach_prepare(u0, p, 0.0)

    f = make_reach_function(prob_f, unflatten_u, unflatten_p, unflatten_t; inplace=true)

    sys = RA.BlackBoxContinuousSystem(f, length(u0)+length(p)+1)
    return new_prob = RA.InitialValueProblem(sys, x0)
end




handle_saturation(f, _) = (f, nothing)
function handle_saturation(f::AbstractSaturated, ::Val{:reachability})
    @unpack f, lb, ub = f
end
