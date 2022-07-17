#=
gauss3.jl
=#


#using ImagePhantoms #: Object, Object3d

export Gauss3
export phantom, radon, spectrum


"""
    Gauss3 <: AbstractShape{3}
"""
struct Gauss3 <: AbstractShape{3} end


# constructors


"""
    Gauss3(cx, cy, cz, wx, wy, wz, Φ=0, Θ=0, value::Number = 1)
    Gauss3(center::NTuple{3,RealU}, radii::NTuple{3,RealU}=(1,1,1), angle::NTuple{2,RealU}=(0,0), v=1)
    Gauss3([9-vector])
    Gauss3(r, v=1) (isotropic of width `w`)
Construct `Gauss3` object from parameters;
here `width` = FWHM (full-width at half-maximum).
"""
function Gauss3(
    cx::RealU,
    cy::RealU,
    cz::RealU,
    wx::RealU,
    wy::RealU = wx,
    wz::RealU = wx,
    Φ::RealU = 0,
    Θ::RealU = 0,
    value::Number = 1,
)
    (cx, cy, cz, wx, wy, wz) = promote(cx, cy, cz, wx, wy, wz)
    Object(Gauss3(), (cx,cy,cz), (wx,wy,wz), (Φ, Θ), value)
end

function Gauss3(
    center::NTuple{3,RealU},
    radii::NTuple{3,RealU} = (1,1,1) .* oneunit(center[1]),
    angle::NTuple{2,RealU} = (0, 0),
    value::Number = 1,
)
    Gauss3(center..., radii..., angle..., value)
end

function Gauss3(v::AbstractVector{<:Number})
    length(v) == 9 || throw(ArgumentError("$v wrong length"))
    Gauss3(v...)
end

Gauss3(w::RealU, v::Number = 1) =
    Gauss3((zero(w),zero(w),zero(w)), (w,w,w), (0,0), v)




# methods


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
    spectrum1(ob::Object3d{Gauss3}, (kx,ky,kz))
Spectrum of unit sphere at `(kx,ky,kz)`, for unitless spatial frequency coordinates.
"""
function spectrum1(ob::Object3d{Gauss3}, kxyz::NTuple{3,Real})
    s = fwhm2spread(1)
    return s^3 * exp(-π * sum(abs2, kxyz) * s^2)
end
