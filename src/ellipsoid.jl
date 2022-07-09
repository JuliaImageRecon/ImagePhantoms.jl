#=
ellipsoid.jl
=#


#using ImagePhantoms #: Object, Object3d

export Ellipsoid
export Sphere
export phantom, radon, spectrum


"""
    Ellipsoid <: AbstractShape3
"""
struct Ellipsoid <: AbstractShape3 end


# constructors


"""
    Ellipsoid(cx, cy, cz, rx, ry, rz, Φ, Θ, value::Number = 1)
    Ellipsoid(center::NTuple{3,RealU}, radii::NTuple{3,RealU}, angle::NTuple{2,RealU}, v)
    Ellipsoid([9-vector])
    Ellipsoid(r, v=1) (sphere of radius `r`)
Construct `Ellipsoid` object from parameters.
"""
function Ellipsoid(
    cx::RealU,
    cy::RealU,
    cz::RealU,
    rx::RealU,
    ry::RealU,
    rz::RealU,
    Φ::RealU = 0,
    Θ::RealU = 0,
    value::Number = 1,
)
    (cx, cy, cz, rx, ry, rz) = promote(cx, cy, cz, rx, ry, rz)
    Object(Ellipsoid(), (cx,cy,cz), (rx,ry,rz), (Φ, Θ), value)
end

function Ellipsoid(
    center::NTuple{3,RealU},
    radii::NTuple{3,RealU},
    angle::NTuple{2,RealU},
    value::Number = 1,
)
    Ellipsoid(center..., radii..., angle..., value)
end

function Ellipsoid(v::AbstractVector{<:Number})
    length(v) == 9 || throw(ArgumentError("$v wrong length"))
    Ellipsoid(v...)
end

Ellipsoid(r::RealU, v::Number = 1) =
    Ellipsoid((zero(r),zero(r),zero(r)), (r,r,r), (0,0), v)


# spheres as a special case

"""
    Sphere(x,y,z,r,v=1) (sphere of radius `r` centered at `(x,y,z)`)
    Sphere((x,y,z), r, v=1) ditto
    Sphere([5-vector]) ditto
    Sphere(r, v=1) centered at origin
Construct `Sphere` objects as special cases of `Ellipsoid` objects.
"""
Sphere(r::RealU, v::Number = 1) = Ellipsoid(r, v)
Sphere(cx::RealU, cy::RealU, cz::RealU, r::RealU, v::Number = 1) =
    Ellipsoid(cx, cy, cz, r, r, r, 0, 0, v)
Sphere(center::NTuple{3,RealU}, r::RealU, v::Number = 1) =
    Ellipsoid(center, (r, r, r), (0, 0), v)

function Sphere(v::AbstractVector{<:Number})
    length(v) == 5 || throw(ArgumentError("$v wrong length"))
    Sphere(v...)
end


# methods


"""
    phantom1(ob::Object3d{Ellipsoid}, (x,y,z))
Evaluate unit sphere at `(x,y,z)`, for unitless coordinates.
"""
phantom1(ob::Object3d{Ellipsoid}, xyz::NTuple{3,Real}) = sum(abs2, xyz) ≤ 1


# x-ray transform (line integral) of unit sphere
# `u,v` should be unitless
function xray1(
    ::Ellipsoid,
    u::Real,
    v::Real,
    ϕ::RealU, # irrelevant
    θ::RealU, # irrelevant
)
    T = promote_type(eltype(u), eltype(v), Float32)
    r2 = u^2 + v^2
    return r2 < 1 ? 2 * sqrt(one(T) - r2) : zero(T)
end


# spectrum

"""
   sphere_transform(f::Real)

Fourier transform of unit-radius sphere.
The argument `f` is the radial coordinate in k-space and is unitless.
See p253 of Bracewell 1978, The Fourier transform and its applications,
or http://doi.org/10.1002/mrm.21292.

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
    spectrum1(ob::Object3d{Ellipsoid}, (kx,ky,kz))
Spectrum of unit sphere at `(kx,ky,kz)`, for unitless spatial frequency coordinates.
"""
function spectrum1(ob::Object3d{Ellipsoid}, kxyz::NTuple{3,Real})
    return sphere_transform(sqrt(sum(abs2, kxyz)))
end
