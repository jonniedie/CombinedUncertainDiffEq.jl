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
    # TODO: Handle integrator, affect!, and affect_neg!
    kwargs = if :callback in keys(prob.kwargs)
        cb = prob.kwargs[:callback]
        condition = (x,t,integrator) -> cb.condition(ComponentArray(x, ax_x), t, integrator)
        @set! cb.condition = condition
        (; prob.kwargs..., callback=cb)
    else
        prob.kwargs
    end
    return ODEProblem(new_f!, getdata(u0), prob.tspan, getdata(p); kwargs...)
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
    new_prob = uncomponentarrayize(prob; u0=u0, p=p)
    return koopman_expectation(g, new_prob, getdata(u0), getdata(p), args...; kwargs...)
end

"""
    mc_solve(prob, u0, p, args...; trajectories, kwargs...)
"""
function mc_solve(prob, u0, p, args...; trajectories, kwargs...)
    u0 = recursive_convert(to_distribution, u0)
    p = recursive_convert(to_distribution, p)
    prob_func = function (prob, i, repeat)
        _u0 = _rand.(u0)
        _p = _rand.(p)
        remake(prob, u0=_u0, p=_p)
    end
    mc_prob = EnsembleProblem(prob, args...; prob_func=prob_func, kwargs...)
    return solve(mc_prob, args...; trajectories=trajectories, kwargs...)
end