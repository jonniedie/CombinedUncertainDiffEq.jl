function soft_step(x; halfwidth=0.1)
    α = 1/sqrt(halfwidth)
    return 1 / (1 + exp(-α*(x+5halfwidth)))
    # if x < -halfwidth
    #     zero(x)
    # elseif x > halfwidth
    #     one(x)
    # else
    #     typeof(x)((sinpi(x/halfwidth/2) + 1)/2)
    # end
end
# soft_step(x, α=1e2) = 1 / (1 + exp(-α*x))

≤ₛ(x, y) = soft_step(y - x)

function soft_clamp(x, lo, hi; kwargs...)
    sgn_lo = soft_step(lo-x; kwargs...)
    sgn_hi = soft_step(x-hi; kwargs...)
    lo_part = lo * sgn_lo
    hi_part = hi * sgn_hi
    mid_part = x * (1-sgn_lo) * (1-sgn_hi)
    return lo_part + mid_part + hi_part
end

function soft_between(x, lo, hi; kwargs...)
    sgn_lo = soft_step(x-lo; kwargs...)
    sgn_hi = soft_step(hi-x; kwargs...)
    return sgn_lo * sgn_hi
end