#=
gauss2.jl
=#

export Gauss2, gauss2


"""
    Gauss2 <: AbstractShape{2}
"""
struct Gauss2 <: AbstractShape{2} end


# constructors


"""
    gauss2(cx, cy, wx, wy=wx, ϕ=0, value::Number=1)
    gauss2(center::NTuple{2,RealU}, width::NTuple{2,RealU}=(1,1), ϕ::RealU=0, v=1)
    gauss2(w, v=1) (isotropic of width `w`)
Construct `Object{Gauss2}` from parameters;
here `width` = FWHM (full-width at half-maximum).

In 1D, the formula is `g(x) = exp(-π ((x - cx) / sx)^2)`
where `sx = fwhm2spread(w) = w * sqrt(π / log(16))`,
which, for `cx=0`,  has 1D FT `G(ν) = sx^2 exp(π (sx νx)^2)`.
"""
gauss2(args... ; kwargs...) = Object(Gauss2(), args...; kwargs...)

gauss2(w::RealU, v::Number = 1) = gauss2((zero(w),zero(w)), (w,w), 0, v)


# helper


"""
    s = fwhm2spread(w)
Convert FWHM `w` to equivalent Gaussian spread `s` for ``\\exp(-π (x/s)^2)``.
`exp(-π (fwhm/2/s)^2) = 1/2` means
`fwhm/2/s) = sqrt(log(2)/π)`
`2s/fwhm = sqrt(π/log(2))`
`2s = fwhm * sqrt(π/log(2))`
`s = fwhm * sqrt(π / log(16))`
"""
@inline fwhm2spread(w) = w * sqrt(π / log(16))


# methods


area1(::Gauss2) = fwhm2spread(1)^2

ℓmax1(::Gauss2) = fwhm2spread(1) # max line integral through a unit gauss2


"""
    phantom1(ob::Object2d{Gauss2}, (x,y))
Evaluate unit gauss2 at `(x,y)`,
for unitless coordinates.
"""
phantom1(ob::Object2d{Gauss2}, xy::NTuple{2,Real}) =
    exp(-π * sum(abs2, xy) / fwhm2spread(1)^2)


# x-ray transform (line integral) of unit gauss2
# `r` should be unitless
function xray1(::Gauss2, r::R, ϕ::RealU) where {R <: Real}
    T = promote_type(R, Float32)
    s = fwhm2spread(1)
    return T(s * exp(-π * abs2(r / s)))
end


"""
    spectrum1(::Object2d{Gauss2}, (kx,ky))
Spectrum of unit gauss2 at `(kx,ky)`,
for unitless spatial frequency coordinates.
"""
function spectrum1(::Object2d{Gauss2}, kxy::NTuple{2,Real})
    s = fwhm2spread(1)
    return s^2 * exp(-π * sum(abs2, kxy) * s^2)
end
