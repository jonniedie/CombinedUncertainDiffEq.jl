"""
    to_distribution(x)

Convert to a `Uniform` distribution from an `Interval`.
"""
to_distribution(x::Interval) = Uniform(x.lo, x.hi)
to_distribution(x::Uniform) = x
to_distribution(x) = x

"""
    to_interval(x)

Convert to a `Interval` from an `Uniform` distribution.
"""
to_interval(x::Interval) = x
to_interval(x::Uniform) = Interval(x.a, x.b)
to_interval(x) = x

"""
    recursive_convert(x)

Recursively convert a variable `x` by applying a conversion function `f`.
Mostly to be used  with `to_distribution` and `to_interval`.
"""
recursive_convert(f, x) = x
recursive_convert(f, x::Union{Interval, Distribution}) = f(x)
recursive_convert(f, x::NamedTuple) = map(x->recursive_convert(f, x), x)
recursive_convert(f, x::AbstractArray) = recursive_convert.(f, x)
