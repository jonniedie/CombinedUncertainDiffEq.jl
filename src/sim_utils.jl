abstract type AbstractSaturated <: Function end

"""
    Saturated(f, lb, ub)

Function wrapper that saturates the output at lower bound `lb` and
upper bound `ub`.
"""
Base.@kwdef @concrete struct Saturated <: AbstractSaturated
    f
    lb
    ub
end

(f::Saturated)(args...; kwargs...) = clamp(f.f(args...; kwargs...), f.lb, f.ub)


Base.@kwdef @concrete struct SoftSaturated <: AbstractSaturated
    f
    lb
    ub
    halfwidth
end

(f::SoftSaturated)(args...; kwargs...) = soft_clamp(f.f(args...; kwargs...), f.lb, f.ub; f.halfwidth)


"""
    WithControls(sys_fun, controls::NamedTuple)

ODE function wrapper that applies control inputs through keyword arguments.

## Example
```julia
function f(x, p, t; u)
    A, B = p.A, p.B
    return A*x + B*u
end

ctrl_f = WithControls(f; u = (x,p,t) -> -p.kp*x)

p = (A=2.5, B=1.5, kp=2.5)  # Parameters
x0 = 1.0                    # Initial Condition
tspan = (0.0, 10.0)         # Time span

prob = ODEProblem(ctrl_f, x0, tspan, p)
sol = solve(prob)
```
"""
@concrete struct WithControls <: AbstractSaturated
    sys_fun
    controls <: NamedTuple
end
WithControls(sys_fun; kwargs...) = WithControls(sys_fun, NamedTuple(kwargs))

(sys::WithControls)(dx, x, p, t) = sys.sys_fun(dx, x, p, t; map(f->maybe_apply(f,x,p,t), sys.controls)...)
(sys::WithControls)(x, p, t) = sys.sys_fun(x, p, t; map(f->maybe_apply(f,x,p,t), sys.controls)...)

DifferentialEquations.isinplace(f::WithControls, args...) = isinplace(f.sys_fun, args...)


"""
    maybe_apply(f, args...)

Maybe apply a function, maybe don't. I don't know. I'm not your boss.
"""
maybe_apply(f::Function, args...) = f(args...)
maybe_apply(val, args...) = val


function saturation_callback(sat_fun::AbstractSaturated)
    @unpack f, lb, ub = sat_fun

    function condition(vars, t, integrator)
        
    end
end

# Needed to make the symbolic stuff work
ComponentArrays.ComponentArray(x::Num) = x

reachable_constraint(f, prob::ODEProblem) = reachable_constraint(f, prob.u0, prob.p, prob.tspan[1])
# reachable_constraint(f, prob::RA.InitialValueProblem) = reachable_constraint(f, prob.x0.center.u, prob.x0.center.p, prob.x0.center.t)
function reachable_constraint(f::Function, u0::ComponentArray, p::ComponentArray, t)
    x0 = make_symbolic_x0(u0, p, t)
    exprs = f(x0.u, x0.p, x0.t)
    return HalfSpace(exprs, x0)
end
function reachable_constraint(fs::AbstractArray, u0::ComponentArray, p::ComponentArray, t)
    x0 = make_symbolic_x0(u0, p, t)
    exprs = map(f->f(x0.u, x0.p, x0.t), fs)
    return HPolyhedron(exprs, x0)
end

function make_symbolic_x0(u0, p, t)
    @variables U[1:length(u0)], P[1:length(p)], T
    u_ca = ComponentArray(collect(U), getaxes(u0))
    p_ca = ComponentArray(collect(P), getaxes(p))
    return ComponentArray(u=u_ca, p=p_ca, t=T)
end

function make_reach_function(f, unflatten_u, unflatten_p, unflatten_t; inplace=true)
    return if inplace
        function (dx, x, _, t)
            U = unflatten_u(x)
            dU = unflatten_u(dx)
            P = unflatten_p(x)
            dP = unflatten_p(dx)
            T = unflatten_t(x)

            f(dU, U, P, T) # Update state
            dP .= zero(x[1])    # Params have no dynamics
            dx[end] = one(x[1]) # Update time variable
            return nothing
        end
    else
        function (x, _, t)
            U = unflatten_u(x)
            P = unflatten_p(x)
            T = unflatten_t(x)
            
            dU = f(U, P, T)
            return [dU; zeros(typeof(x[1]), length(P)); one(x[1])]
        end
    end
end

function reach_prepare(u0, p, t)
    nu = length(u0)
    np = length(p)
    u0, unflatten_u = _reach_prepare(u0)
    p, unflatten_p = _reach_prepare(p, nu)
    t, unflatten_t = _reach_prepare(t, nu+np)
    x0 = prod(u0) × prod(p) × t |> concretize
    # center = ComponentArray(u=unflatten_u(x0.center), p=unflatten_p(x0.center), t=unflatten_t(x0.center))
    # radius = ComponentArray(u=unflatten_u(x0.radius), p=unflatten_p(x0.radius), t=unflatten_t(x0.radius))
    # x0 = Hyperrectangle(center, radius)
    return x0, unflatten_u, unflatten_p, unflatten_t
end

function reach_prepare(u0, p)
    nu = length(u0)
    u0, unflatten_u = _reach_prepare(u0)
    p, unflatten_p = _reach_prepare(p, nu)
    x0 = prod(u0) × prod(p) |> concretize
    # center = ComponentArray(u=unflatten_u(x0.center), p=unflatten_p(x0.center))
    # radius = ComponentArray(u=unflatten_u(x0.radius), p=unflatten_p(x0.radius))
    # x0 = Hyperrectangle(center, radius)
    return x0, unflatten_u, unflatten_p
end

function _reach_prepare(p, offset=0)
    p = to_reach_interval(p)
    i = offset + 1
    unflatten = x -> x[i]
    return p, unflatten
end
function _reach_prepare(p::AbstractArray, offset=0)
    np = length(p)
    p = recursive_convert(to_reach_interval, p)
    i = (1:np) .+ offset
    unflatten = x -> view(x, i)
    return p, unflatten
end
function _reach_prepare(p::ComponentArray, offset=0)
    ax_p = getaxes(p)
    np = length(p)
    p = recursive_convert(to_reach_interval, p)
    i = (1:np) .+ offset
    unflatten = x -> ComponentArray(view(x, i), ax_p)
    return p, unflatten
end