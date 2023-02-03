#=
jinc.jl
=#

export jinc

import SpecialFunctions


"""
    jinc(x)

Return `jinc(x) = J1(π*x)/(2x)`,
where `J1` is a Bessel function of the first kind.

The argument `x` must be unitless.

!!! note
    `SpecialFunctions.jinc(0) = 1` whereas the convention here is `jinc(0) = π / 4`
    which corresponds to the area of a disk of unit diameter.
"""
function jinc(x::X) where {X <: Real}
    X(π) / 4 * SpecialFunctions.jinc(x)
end
