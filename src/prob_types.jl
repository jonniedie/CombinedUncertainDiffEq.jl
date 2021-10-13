abstract type AbstractUncertainODEProblem end


@concrete struct MonteCarloProblem <: AbstractUncertainODEProblem
    prob <: EnsembleProblem
    kwargs

    function MonteCarloProblem(prob; u0=prob.u0, p=prob.p, kwargs...)
        u0 = recursive_convert(to_distribution, u0)
        p = recursive_convert(to_distribution, p)
        prob_func = function (prob, i, repeat)
            _u0 = _rand.(u0)
            _p = _rand.(p)
            return remake(prob, u0=_u0, p=_p)
        end
        ens_prob = EnsembleProblem(prob; prob_func=prob_func)
        return new{typeof(ens_prob), typeof(kwargs)}(ens_prob, kwargs)
    end
end


@concrete struct KoopmanProblem <: AbstractUncertainODEProblem
    prob <: ODEProblem
    u0_CoV <: Function
    p_CoV <: Function
    kwargs

    function KoopmanProblem(prob::ODEProblem; u0=prob.u0, p=prob.p, kwargs...)
        u0_CoV, p_CoV = make_conversion_functions(u, p)
        u0 = recursive_convert(to_distribution, getdata(u0))
        p = recursive_convert(to_distribution, getdata(p))
        prob = remake(prob; u0=u0, p=p)
        return new{typeof(prob), typeof(u0_CoV), typeof(p_CoV)}(prob, u0_CoV, p_CoV, kwargs)
    end
end
KoopmanProblem(prob::KoopmanProblem; u0=prob.u0, p=prob.p, kwargs...) = KoopmanProblem(remake(prob.prob; u0=u0, p=p); prob.kwargs..., kwargs...)
KoopmanProblem(f::Function, args...; kwargs...) = KoopmanProblem(ODEProblem(f, args...; kwargs...))


make_conversion_functions(u0, p) = ((u,p)->u, (u,p)->p)
function make_conversion_functions(u0::ComponentArray, p::ComponentArray)
    ax_u0 = getaxes(u0)
    ax_p = getaxes(p)
    u0_CoV = (u, p) -> ComponentArray(u, ax_u0)
    p_CoV = (u, p) -> ComponentArray(p, ax_p)
    return u0_CoV, p_CoV
end