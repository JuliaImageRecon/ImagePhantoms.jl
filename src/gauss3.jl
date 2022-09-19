#=
gauss3.jl
=#

export Gauss3, gauss3


"""
    Gauss3 <: AbstractShape{3}
"""
struct Gauss3 <: AbstractShape{3} end


# constructor


"""
    gauss3(cx, cy, cz, wx, wy, wz, Φ=0, Θ=0, value::Number = 1)
    gauss3(center::NTuple{3,RealU}, radii::NTuple{3,RealU}=(1,1,1), angle::NTuple{2,RealU}=(0,0), v=1)
    gauss3(r, v=1) (isotropic of width `w`)
Construct `Object{Gauss3}` from parameters;
here `width` = FWHM (full-width at half-maximum).
"""
gauss3(args... ; kwargs...) = Object(Gauss3(), args...; kwargs...)


# methods


volume1(::Gauss3) = fwhm2spread(1)^3

ℓmax1(::Gauss3) = fwhm2spread(1) # max line integral through a unit gauss2


"""
    phantom1(ob::Object3d{Gauss3}, (x,y,z))
Evaluate Gauss3 `(x,y,z)`, for unitless coordinates.
"""
phantom1(ob::Object3d{Gauss3}, xyz::NTuple{3,Real}) =
        exp(-π * sum(abs2, xyz) / fwhm2spread(1)^2)


# x-ray transform (line integral) of "unit" Gauss3
# `u,v` should be unitless
function xray1(
    ::Gauss3,
    u::Real,
    v::Real,
    ϕ::RealU, # irrelevant
    θ::RealU, # irrelevant
)
    s = fwhm2spread(1)
    return s * exp(-π * (u^2 + v^2) / s^2)
end


"""
    spectrum1(::Object3d{Gauss3}, (kx,ky,kz))
Spectrum of unit sphere at `(kx,ky,kz)`, for unitless spatial frequency coordinates.
"""
function spectrum1(::Object3d{Gauss3}, kxyz::NTuple{3,Real})
    s = fwhm2spread(1)
    return s^3 * exp(-π * sum(abs2, kxyz) * s^2)
end
