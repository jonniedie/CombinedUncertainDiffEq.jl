"""
    to_distribution(x)

Convert to a `Uniform` distribution from an `Interval`.
"""
to_distribution(x::Interval) = Uniform(x.lo, x.hi)
to_distribution(x::Distribution) = x
to_distribution(x::ReachInterval) = to_distribution(x.dat)
to_distribution(x) = x


"""
    to_interval(x)

Convert to a `Interval` from a `Distribution`.
"""
to_interval(x::Interval) = x
to_interval(x::Distribution) = Interval(_minimum(x), _maximum(x))
to_interval(x::ReachInterval) = x.dat
to_interval(x) = x


"""
    to_reach_interval(x)

Convert to a `ReachInterval` from an `Interval` or `Distribution`.
"""
to_reach_interval(x::Interval) = ReachInterval(x)
to_reach_interval(x::Distribution) = to_reach_interval(to_interval(x))
to_reach_interval(x::ReachInterval) = x
to_reach_interval(x) = RA.Singleton([x])
to_reach_interval(x::Integer) = to_reach_interval(float(x))


"""
    recursive_convert(x)

Recursively convert a variable `x` by applying a conversion function `f`.
Mostly to be used  with `to_distribution` and `to_interval`.
"""
recursive_convert(f, x) = f(x)
recursive_convert(f, x::NamedTuple) = map(x->recursive_convert(f, x), x)
recursive_convert(f, x::AbstractArray) = recursive_convert.(f, x)


# Internal maximum and minimum to work across all uncertainty types
_minimum(x) = minimum(x)
_minimum(x::Interval) = inf(x)
_minimum(x::ReachInterval) = min(x)

_maximum(x) = maximum(x)
_maximum(x::Interval) = sup(x)
_maximum(x::ReachInterval) = max(x)