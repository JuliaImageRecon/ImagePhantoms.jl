#=
ellipsoid.jl
=#

export Ellipsoid, ellipsoid, sphere


"""
    Ellipsoid <: AbstractShape{3}
"""
struct Ellipsoid <: AbstractShape{3} end


# constructor


"""
    ellipsoid(cx, cy, cz, rx=1, ry=1, rz=1, Φ=0, Θ=0, value::Number = 1)
    ellipsoid(center::NTuple{3,RealU}, radii::NTuple{3,RealU}, angle::NTuple{3,RealU}, v)
Construct `Object{Ellipsoid}` from parameters.
"""
ellipsoid(args... ; kwargs...) = Object(Ellipsoid(), args...; kwargs...)


# special case: spheres

"""
    sphere(x, y, z, r,v=1) (sphere of radius `r` centered at `(x,y,z)`)
    sphere((x,y,z), r, v=1) ditto
    sphere(r, v=1) centered at origin
Construct spheres as special cases of `Ellipsoid`.
"""
sphere(cx::RealU, cy::RealU, cz::RealU, r::RealU, v::Number = 1) =
    ellipsoid(cx, cy, cz, r, r, r, 0, 0, v)
sphere(center::NTuple{3,RealU}, r::RealU, v::Number = 1) =
    sphere(center..., r, v)
sphere(r::RealU, v::Number = 1) = sphere((zero(r), zero(r), zero(r)), v)


# methods


volume1(::Ellipsoid) = 4/3 * π # volume of unit sphere

ℓmax1(::Ellipsoid) = 2 # maximum chord through a unit sphere


"""
    phantom1(ob::Object3d{Ellipsoid}, (x,y,z))
Evaluate unit sphere at `(x,y,z)`,
for unitless coordinates.
"""
phantom1(ob::Object3d{Ellipsoid}, xyz::NTuple{3,Real}) = (sum(abs2, xyz) ≤ 1)


# x-ray transform (line integral) of unit sphere
# `u,v` should be unitless
function xray1(
    ::Ellipsoid,
    u::Ru,
    v::Rv,
    ϕ::RealU, # irrelevant
    θ::RealU, # irrelevant
) where {Ru <: Real, Rv <: Real}
    T = promote_type(Ru, Rv, Float32)
    r2 = u^2 + v^2
    return (r2 < 1 ? 2 * sqrt(one(T) - r2) : zero(T))::T
end


# spectrum

"""
   sphere_transform(f::Real)

Fourier transform of unit-radius sphere.
The argument `f` is the radial coordinate in k-space and is unitless.
See p253 of Bracewell 1978, The Fourier transform and its applications,
or https://doi.org/10.1002/mrm.21292.

Formula: `4/3 π` for `f ≈ 0`, otherwise
`(sin(2πf) - 2πf cos(2πf)) / (2 * π^2 * f^3)`.
"""
function sphere_transform(f::T) where {T <: AbstractFloat}
    atol = eps(T)
    if abs(f)^3 ≤ atol
         return 4/3 * π
    end
    f2pi = 2π * f
    (s, c) = sincos(f2pi)
    # numerator: s - f2pi * c ≈ (2πf)^3/3 for small f
    return (s - f2pi * c) / (2 * π^2 * f^3)
end


"""
    spectrum1(::Object3d{Ellipsoid}, (kx,ky,kz))
Spectrum of unit sphere at `(kx,ky,kz)`,
for unitless spatial frequency coordinates.
"""
function spectrum1(::Object3d{Ellipsoid}, kxyz::NTuple{3,Real})
    return sphere_transform(sqrt(sum(abs2, kxyz)))
end
