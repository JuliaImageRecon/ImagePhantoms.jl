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


# helper

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




# methods


"""
    phantom1(ob::Object3d{Ellipsoid}, (x,y,z))
Evaluate unit sphere at `(x,y,z)`, for unitless coordinates.
"""
phantom1(ob::Object3d{Ellipsoid}, xyz::NTuple{3,Real}) = sum(abs2, xyz) ≤ 1


#=
"""
    xray_ellipsoid(s, t, ϕ, θ, cx, cy, cz, rx, ry, rz, Φ, Θ)
X-ray transform at `(s, t, ϕ, θ)` of ellipsoid.
"""
function xray_ellipsoid(s, t, ϕ, θ, cx, cy, rx, ry, rz, Φ, Θ)
    Θ == 0 || throw("polar rotation not done")
    (sinϕ, cosϕ) = sincos(ϕ)
    r -= cx * cosϕ + cy * sinϕ # Radon translation property
    (sinϕ, cosϕ) = sincos(ϕ - θ) # Radon rotation property
    return lmax * trapezoid(r, -dmax, -dbreak, dbreak, dmax)
end
=#


"""
Translated from ellipsoid_proj.m in MIRT
Caution: note that `ϕ, θ` denote projection view angles
whereas `Φ, Θ` denote object rotation angles.
"""
function xray_ellipsoid(u, v, ϕ, polar, cx, cy, cz, rx, ry, rz, xang, zang)
#function xray_ellipsoid(u, v, β, θ, cx, cy, cz, rx, ry, rz, xang, zang)
#function xray_ellipsoid(u, v, β, θ, cx, cy, cz, rx, ry, rz, Φ, Θ)
#function ellipsoid_proj_do(ss, tt, beta, source_zs, dso, dod, dfs, oversample)

    zang == 0 || throw("Z angle not done")
    (spolar, cpolar) = sincos(polar)

    sinaz, cosaz = sincos(ϕ)
    ushift = cx * cosaz + cy * sinaz
    vshift = (cx * sinaz - cy * cosaz) * spolar + cz * cpolar
    u -= ushift
    v -= vshift

    az = ϕ - xang
    sinaz, cosaz = sincos(az)
    p1 = u * cosaz + v * sinaz * spolar
    p2 = u * sinaz - v * cosaz * spolar
    p3 = v * cpolar

    e1 = -sinaz * cpolar
    e2 = cosaz * cpolar
    e3 = spolar

    A = (e1 / rx)^2 + (e2 / ry)^2 + (e3 / rz)^2
    B = p1 * e1 / rx^2 + p2 * e2 / ry^2 + p3 * e3 / rz^2
    C = (p1 / rx)^2 + (p2 / ry)^2 + (p3 / rz)^2 - 1

    return 2 * real(sqrt(Complex(B^2 - A*C))) / A
end


"""
    radon(ob::Object3d{Ellipsoid})
Returns function of `(s,t,ϕ,θ)` for evaluating a line integral of an Ellipsoid,
where `ϕ` is the azimuthal angle and `θ` is the polar angle.
"""
radon(ob::Object3d{Ellipsoid}) = (s,t,ϕ,θ) -> ob.value *
    xray_ellipsoid(s, t, ϕ, θ, ob.center..., ob.width..., ob.angle...)


function spectrum_ellipsoid(fx, fy, fz, cx, cy, cz, rx, ry, rz, Φ, Θ)
    (kx, ky, kz) = rotate3d(fx, fy, fz, Φ, Θ) # rotate then translate
    return rx * ry * rz * cispi(-2*(fx*cx + fy*cy + fz*cz)) *
        sphere_transform(sqrt((rx*kx)^2 + (ry*ky)^2 + (rz*kz)^2))
end


"""
    spectrum(ob::Object3d{Ellipsoid})
Returns function of ``(f_x,f_y,f_z)`` for the spectrum (3D Fourier transform).
"""
spectrum(ob::Object3d{Ellipsoid}) = (fx,fy,fz) -> ob.value *
    spectrum_ellipsoid(fx, fy, fz, ob.center..., ob.width..., ob.angle...)
