#=
cuboid.jl
actually a rectangular cuboid https://en.wikipedia.org/wiki/Cuboid
=#


#using ImagePhantoms #: Object, Object3d

export Cuboid, cuboid, cube
export phantom, radon, spectrum


"""
    Cuboid <: AbstractShape{3}
"""
struct Cuboid <: AbstractShape{3} end


# constructor


"""
    cuboid(cx, cy, cz, wx, wy, wz, Φ, Θ, value::Number)
    cuboid(center::NTuple{3,RealU}, width::NTuple{3,RealU}, angle::NTuple{2,RealU}, v)
    cuboid([9-vector])
Construct `Object{Cuboid}` from parameters;
here `width` is the full-width.
"""
cuboid(args... ; kwargs...) = Object(Cuboid(), args...; kwargs...)


# special case: cubes

"""
    cube(x, y, z, w, v=1) (cube of width `w` centered at `(x,y,z)`)
    cube((x,y,z), w, v=1) ditto
    cube([5-vector]) ditto
    cube(w, v=1) centered at origin
Construct cubes as special cases of `Cuboid`.
"""
cube(cx::RealU, cy::RealU, cz::RealU, w::RealU, v::Number = 1) =
    cuboid(cx, cy, cz, w, w, w, 0, 0, v)
cube(center::NTuple{3,RealU}, w::RealU, v::Number = 1) =
    cube(center..., w, v)
cube(w::RealU, v::Number = 1) = cube((zero(w), zero(w), zero(w)), v)

function cube(v::AbstractVector{<:Number})
    length(v) == 5 || throw(ArgumentError("$v wrong length"))
    cube(v...)
end


# methods


"""
    phantom1(ob::Object3d{Cuboid}, (x,y,z))
Evaluate unit cube at `(x,y,z)`,
for unitless coordinates.
"""
phantom1(ob::Object3d{Cuboid}, xyz::NTuple{3,Real}) = (maximum(abs, xyz) ≤ 0.5)


# radon


_sort(a, b) = a < b ? (a, b) : (b, a)

"""
    (ℓmin, ℓmax) = cube_bounds(p, e)
Bounds of ℓ corresponding to rect(x) for direction `e` from point `p`.
"""
function cube_bounds(p::T, e::T) where T <: AbstractFloat
    return abs(e) < eps(T) ? (-T(Inf), T(Inf)) :
        _sort((-T(0.5) - p) / e, ( T(0.5) - p) / e)
end
function cube_bounds(p::Real, e::Real)
   T = promote_type(eltype(p), eltype(e), Float32)
   return cube_bounds(T(p), T(e))
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

    e1 = -sϕ * cθ # x = p1 + ℓ * e1
    e2 = cϕ * cθ  # y = p2 + ℓ * e2
    e3 = sθ       # z = p3 + ℓ * e3

    ℓxmin, ℓxmax = cube_bounds(p1, e1)
    ℓymin, ℓymax = cube_bounds(p2, e2)
    ℓzmin, ℓzmax = cube_bounds(p3, e3)

    minℓ = max(ℓxmin, ℓymin, ℓzmin)
    maxℓ = min(ℓxmax, ℓymax, ℓzmax)
    ℓ = max(maxℓ - minℓ, zero(T))

    if abs(e1) < eps(T)
        ℓ *= (-1/2 ≤ u < 1/2)
    end
    if abs(e2) < eps(T)
        ℓ *= (-1/2 ≤ u < 1/2)
    end
    if abs(e3) < eps(T)
        ℓ *= (-1/2 ≤ v < 1/2)
    end
    return ℓ
end


# spectrum

"""
    spectrum(ob::Object3d{Cuboid}, (kx,ky,kz))
Spectrum of unit cube at `(kx,ky,kz)`,
for unitless spatial frequency coordinates.
"""
function spectrum1(ob::Object3d{Cuboid}, kxyz::NTuple{3,Real})
    return prod(sinc, kxyz)
end
