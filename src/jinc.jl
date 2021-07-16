#=
jinc.jl
=#

export jinc

using SpecialFunctions: besselj1


"""
    jinc(x)

Return `jinc(x) = J1(π*x)/(2x)`,
where `J1` is a Bessel function of the first kind.

The argument `x` must be unitless.

Return type is `promote_type(typeof(x), Float32)`.
"""
function jinc(x::X) where {X <: Real}
    T = promote_type(X, Float32)
    if (x == 0) return convert(T, π/4) end
    y = abs(x)
    convert(T, besselj1(π*y) / (2*y))
end
