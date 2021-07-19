#=
cuboid.jl
actually a rectangular cuboid https://en.wikipedia.org/wiki/Cuboid
=#


#using ImagePhantoms #: Object, Object3d

export Cuboid
export Cube
export phantom, radon, spectrum


"""
    Cuboid <: AbstractShape3
"""
struct Cuboid <: AbstractShape3 end


# constructors


"""
    Cuboid(cx, cy, cz, wx, wy, wz, Φ, Θ, value::Number)
    Cuboid(center::NTuple{3,RealU}, width::NTuple{3,RealU}, angle::NTuple{2,RealU}, v)
    Cuboid([9-vector])
    Cuboid(r, v=1) (cube of width `r`)
Construct `Cuboid` object from parameters;
here `width` is the full-width.
"""
function Cuboid(
    cx::RealU,
    cy::RealU,
    cz::RealU,
    wx::RealU,
    wy::RealU,
    wz::RealU,
    Φ::RealU = 0,
    Θ::RealU = 0,
    value::Number = 1,
)
    (cx, cy, cz, wx, wy, wz) = promote(cx, cy, cz, wx, wy, wz)
    Object(Cuboid(), (cx,cy,cz), (wx,wy,wz), (Φ, Θ), value)
end

function Cuboid(
    center::NTuple{3,RealU},
    width::NTuple{3,RealU},
    angle::NTuple{2,RealU},
    value::Number = 1,
)
    Cuboid(center..., width..., angle..., value)
end

function Cuboid(v::AbstractVector{<:Number})
    length(v) == 9 || throw(ArgumentError("$v wrong length"))
    Cuboid(v...)
end

Cuboid(r::RealU, v::Number = 1) =
    Cuboid((zero(r),zero(r),zero(r)), (r,r,r), (0,0), v)


# cubes as a special case

"""
    Cube(x,y,z,w,v=1) (cube of width `w` centered at `(x,y,z)`)
    Cube((x,y,z), w, v=1) ditto
    Cube([5-vector]) ditto
    Cube(w, v=1) centered at origin
Construct `Cube` objects as special cases of `Cuboid` objects.
"""
Cube(w::RealU, v::Number = 1) = Cuboid(w, v)
Cube(cx::RealU, cy::RealU, cz::RealU, w::RealU, v::Number = 1) =
    Cuboid(cx, cy, cz, w, w, w, 0, 0, v)
Cube(center::NTuple{3,RealU}, w::RealU, v::Number = 1) =
    Cuboid(center, (w, w, w), (0, 0), v)

function Cube(v::AbstractVector{<:Number})
    length(v) == 5 || throw(ArgumentError("$v wrong length"))
    Cube(v...)
end


# helper

#=
"""
    (lxmin, lxmax) = cuboid_proj_line1(rx, p1, e1)
"""
function cuboid_proj_line1(rx, p1, e1)
    tmp = (e1 == 0)
    e1(tmp) = inf
    # bounds of l corresponding to rect(x/rx)
    lxmin = (-rx/2 - p1) ./ e1
    lxmax = ( rx/2 - p1) ./ e1
    # re-arrange the bounds so that lxmin contains the minimum l values
    # and lxmax contains the maximum l values
    temp = lxmin
    lxmin = min(lxmin, lxmax)
    lxmax = max(temp, lxmax)
    # exclude points where e1=0 by setting lxmin = -Inf and lxmax = Inf
    lxmin(tmp) = -inf
    lxmax(tmp) = inf
    return (lxmin, lxmax)
end
=#


# methods


"""
    phantom(ob::Object3d{Cuboid})
Returns function of `(x,y)` for making image.
"""
phantom(ob::Object3d{Cuboid}) = (x,y,z) ->
    ob.value * (maximum(abs, coords(ob, x, y, z)) ≤ 0.5)


#=
todo, use cuboid_proj.m
"""
    xray_cuboid(s, t, ϕ, θ, cx, cy, cz, wx, wy, wz, Φ, Θ)
X-ray transform at `(s, t, ϕ, θ)` of cuboid.
"""
function xray_cuboid(s, t, ϕ, θ, cx, cy, wx, wy, Φ, Θ)
    Θ == 0 || throw("polar rotation not done")
    (sinϕ, cosϕ) = sincos(ϕ)
    r -= cx * cosϕ + cy * sinϕ # Radon translation property
    (sinϕ, cosϕ) = sincos(ϕ - θ) # Radon rotation property
    xmax = wx * abs(cosϕ)
    ymax = wy * abs(sinϕ)
    lmax = wx * wy / max(xmax, ymax)
    dmax = (xmax + ymax) / 2
    dbreak = abs(xmax - ymax) / 2
    return lmax * trapezoid(r, -dmax, -dbreak, dbreak, dmax)
end
=#

"""
    radon(ob::Object3d{Cuboid})
Returns function of `(s,t,ϕ,θ)` for evaluating the projection of a Cuboid.
"""
radon(ob::Object3d{Cuboid}) = (s,t,ϕ,θ) -> ob.value *
    xray_cuboid(s, t, ϕ, θ, ob.center..., ob.width..., ob.angle...)


function spectrum_cuboid(fx, fy, fz, cx, cy, cz, wx, wy, Φ, Θ)
    (kx, ky, kz) = rotate3d(fx, fy, fz, Φ, Θ) # rotated, then translate
    return wx * sinc(kx * wx) * exp(-2im*π*fx*cx) *
           wy * sinc(ky * wy) * exp(-2im*π*fy*cy) *
           wz * sinc(kz * wz) * exp(-2im*π*fz*cz)
end

"""
    spectrum(ob::Object3d{Cuboid})
Returns function of ``(f_x,f_y,f_z)`` for the spectrum (3D Fourier transform)
of a Cuboid.
"""
spectrum(ob::Object3d{Cuboid}) = (fx,fy,fz) -> ob.value *
    spectrum_cuboid(fx, fy, fz, ob.center..., ob.width..., ob.angle...)
