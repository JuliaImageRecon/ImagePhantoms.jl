#=
gauss2.jl
=#


#using ImagePhantoms #: Object, Object2d

export Gauss2
export phantom, radon, spectrum, fwhm2spread


"""
    Gauss2 <: AbstractShape2
"""
struct Gauss2 <: AbstractShape2 end


# constructors


"""
    Gauss2(cx, cy, wx, wy=wx, ϕ=0, value::Number=1)
    Gauss2(center::NTuple{2,RealU}, width::NTuple{2,RealU}=(1,1), ϕ::RealU=0, v=1)
    Gauss2([6-vector])
    Gauss2(w, v=1) (isotropic of width `w`)
Construct `Gauss2` object from parameters;
here `width` = FWHM (full-width at half-maximum).

In 1D, the formula is `g(x) = exp(-π ((x - cx) / sx)^2)`
where `sx = fwhm2spread(w) = w * sqrt(π / log(16))`,
which, for `cx=0`,  has 1D FT `G(ν) = sx^2 exp(π (sx νx)^2)`.
"""
function Gauss2(
    cx::RealU,
    cy::RealU,
    wx::RealU,
    wy::RealU = wx,
    ϕ::RealU = 0,
    value::Number = 1,
)
    (cx, cy, wx, wy) = promote(cx, cy, wx, wy)
    Object(Gauss2(), (cx,cy), (wx,wy), ϕ, value)
end

function Gauss2(
    center::NTuple{2,RealU},
    width::NTuple{2,RealU} = (1,1) .* oneunit(center[1]),
    ϕ::RealU = 0,
    value::Number = 1,
)
    Gauss2(center..., width..., ϕ, value)
end

function Gauss2(v::AbstractVector{<:Number})
    length(v) == 6 || throw(ArgumentError("$v wrong length"))
    Gauss2(v...)
end

Gauss2(w::RealU, v::Number = 1) = Gauss2((zero(w),zero(w)), (w,w), 0, v)


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


"""
    phantom1(ob::Object2d{Gauss2}, (x,y))
Evaluate unit gauss2 at `(x,y)`,
for unitless coordinates.
"""
phantom1(ob::Object2d{Gauss2}, xy::NTuple{2,Real}) =
        exp(-π * sum(abs2, xy) / fwhm2spread(1)^2)


# x-ray transform (line integral) of unit gauss2
# `r` should be unitless
function xray1(::Gauss2, r::Real, ϕ::RealU)
    s = fwhm2spread(1)
    return s * exp(-π * abs2(r / s))
end


"""
    spectrum1(ob::Object2d{Gauss2}, (kx,ky))
Spectrum of unit gauss2 at `(kx,ky)`,
for unitless spatial frequency coordinates.
"""
function spectrum1(ob::Object2d{Gauss2}, kxy::NTuple{2,Real})
    s = fwhm2spread(1)
    return s^2 * exp(-π * sum(abs2, kxy) * s^2)
end
