#=
gauss2.jl
=#


#using ImagePhantoms #: Object, Object2d

export Gauss2
export phantom, radon, spectrum


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
    σ = fwhm2sigma(w)
Convert FWHM `w` to equivalent Gaussian `σ` for ``\\exp(-π (x/σ)^2)``.
"""
@inline fwhm2sigma(w) = w / sqrt(log(256))


# methods


"""
    phantom(ob::Object2d{Gauss2})
Returns function of `(x,y)` for making image.
"""
function phantom(ob::Object2d{Gauss2})
    return (x,y) -> ob.value * # trick due to fwhm to σ scaling:
        exp(-π * sum(abs2.(coords(ob, x, y))) / fwhm2sigma(1)^2)
end


"""
    radon_gauss2(r, ϕ, cx, cy, wx, wy, θ)
Radon transform at `(r,ϕ)` of 2D Gaussian.
"""
function radon_gauss2(r, ϕ, cx, cy, wx, wy, θ)
    (σx, σy) = fwhm2sigma.((wx, wy))
    (sinϕ, cosϕ) = sincos(ϕ)
    r -= cx * cosϕ + cy * sinϕ # Radon translation property
    (sinϕ, cosϕ) = sincos(ϕ - θ) # Radon rotation property
    σ = sqrt(abs2(σx * cosϕ) + abs2(σy * sinϕ)) # by Fourier-slice Thm.
    return σx * σy / σ * exp(-π * abs2(r / σ))
end


"""
    radon(ob::Object2d{Gauss2})
Returns function of `(r,ϕ)` for making a sinogram.
"""
radon(ob::Object2d{Gauss2}) = (r,ϕ) -> ob.value *
    radon_gauss2(r, ϕ, ob.center..., ob.width..., ob.angle[1])


function spectrum_gauss2(fx, fy, cx, cy, wx, wy, θ)
    (σx, σy) = fwhm2sigma.((wx, wy))
    (kx, ky) = rotate2d(fx, fy, θ) # rotate first, then translate
    return σx * exp(-2im*π*fx*cx) *
           σy * exp(-2im*π*fy*cy) * exp(-π * (abs2(σx*kx) + abs2(σy*ky)))
end

"""
    spectrum(ob::Object2d{Gauss2})
Returns function of ``(f_x,f_y)`` for the spectrum (2D Fourier transform).
"""
spectrum(ob::Object2d{Gauss2}) = (fx,fy) -> ob.value *
    spectrum_gauss2(fx, fy, ob.center..., ob.width..., ob.angle[1])
