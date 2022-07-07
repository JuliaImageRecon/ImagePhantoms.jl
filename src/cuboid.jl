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


# methods


"""
    phantom1(ob::Object3d{Cuboid}, (x,y,z))
Evaluate unit cube at `(x,y,z)`, for unitless coordinates.
"""
phantom1(ob::Object3d{Cuboid}, xyz::NTuple{3,Real}) = (maximum(abs, xyz) ≤ 0.5)


# radon

"""
    (ℓmin, ℓmax) = cuboid_proj_line1(p1, e1)
"""
function cuboid_proj_line1(p1, e1)
    if e1 == 0
        return (-Inf,Inf)
    end
    # bounds of ℓ corresponding to rect(x)
    ℓmin = (-1/2 - p1) ./ e1
    ℓmax = ( 1/2 - p1) ./ e1
    return min(ℓmin, ℓmax), max(ℓmin, ℓmax)
end


# x-ray transform (line integral) of unit cube
# `u,v` should be unitless
function xray1(
    ::Cuboid,
    u::Real,
    v::Real,
    ϕ::RealU, # azim
    θ::RealU, # polar
)
    T = promote_type(eltype(u), eltype(v), Float32)

    (sϕ, cϕ) = sincos(ϕ)
    (sθ, cθ) = sincos(θ)

    p1 = u * cϕ + v * sϕ * sθ
    p2 = u * sϕ - v * cϕ * sθ
    p3 = v * cθ

    e1 = -sϕ * cθ # x = p1 + l*e1
    e2 = cϕ * cθ  # y = p2 + l*e2
    e3 = sθ       # z = p3 + l*e3

    ℓxmin, ℓxmax = cuboid_proj_line1(p1, e1)
    ℓymin, ℓymax = cuboid_proj_line1(p2, e2)
    ℓzmin, ℓzmax = cuboid_proj_line1(p3, e3)

    ℓmin = max(ℓxmin, ℓymin, ℓzmin)
    ℓmax = min(ℓxmax, ℓymax, ℓzmax)
    ℓ = max(ℓmax - ℓxmin, zero(T))
    if e1 == 0 && !(-1/2 ≤ u ≤ 1/2)
        return zero(T)
    end
    if e2 == 0 && !(-1/2 ≤ u ≤ 1/2)
        return zero(T)
    end
    if e3 == 0 && !(-1/2 ≤ v ≤ 1/2)
        return zero(T)
    end
    return ℓ
end


# spectrum

"""
    spectrum(ob::Object3d{Cuboid}, (kx,ky,kz))
Spectrum of unit cube at `(kx,ky,kz)`, for unitless spatial frequency coordinates.
"""
function spectrum1(ob::Object3d{Cuboid}, kxyz::NTuple{3,Real})
    return prod(sinc, kxyz)
end
