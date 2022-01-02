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
    phantom(ob::Object2d{Gauss2})
Returns function of `(x,y)` for making image.
"""
function phantom(ob::Object2d{Gauss2})
    return (x,y) -> ob.value * # trick due to fwhm to spread scaling:
        exp(-π * sum(abs2.(coords(ob, x, y))) / fwhm2spread(1)^2)
end


"""
    radon_gauss2(r, ϕ, cx, cy, wx, wy, θ)
Radon transform at `(r,ϕ)` of 2D Gaussian.
"""
function radon_gauss2(r, ϕ, cx, cy, wx, wy, θ)
    (sx, sy) = fwhm2spread.((wx, wy))
    (sinϕ, cosϕ) = sincos(ϕ)
    r -= cx * cosϕ + cy * sinϕ # Radon translation property
    (sinϕ, cosϕ) = sincos(ϕ - θ) # Radon rotation property
    s = sqrt(abs2(sx * cosϕ) + abs2(sy * sinϕ)) # by Fourier-slice Thm.
    return sx * sy / s * exp(-π * abs2(r / s))
end


"""
    radon(ob::Object2d{Gauss2})
Returns function of `(r,ϕ)` for making a sinogram.
"""
radon(ob::Object2d{Gauss2}) = (r,ϕ) -> ob.value *
    radon_gauss2(r, ϕ, ob.center..., ob.width..., ob.angle[1])


function spectrum_gauss2(fx, fy, cx, cy, wx, wy, θ)
    (sx, sy) = fwhm2spread.((wx, wy))
    (kx, ky) = rotate2d(fx, fy, θ) # rotate first, then translate
    return sx * exp(-2im*π*fx*cx) *
           sy * exp(-2im*π*fy*cy) * exp(-π * (abs2(sx*kx) + abs2(sy*ky)))
end

"""
    spectrum(ob::Object2d{Gauss2})
Returns function of ``(f_x,f_y)`` for the spectrum (2D Fourier transform).
"""
spectrum(ob::Object2d{Gauss2}) = (fx,fy) -> ob.value *
    spectrum_gauss2(fx, fy, ob.center..., ob.width..., ob.angle[1])
