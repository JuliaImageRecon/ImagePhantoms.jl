#=
cuboid.jl
actually a rectangular cuboid https://en.wikipedia.org/wiki/Cuboid
=#


export Cuboid, cuboid, cube


"""
    Cuboid <: AbstractShape{3}
"""
struct Cuboid <: AbstractShape{3} end


# constructor


"""
    cuboid(cx, cy, cz, wx, wy, wz, Φ, Θ, value::Number)
    cuboid(center::NTuple{3,RealU}, width::NTuple{3,RealU}, angle::NTuple{3,RealU}, v)
Construct `Object{Cuboid}` from parameters;
here `width` is the full-width.
"""
cuboid(args... ; kwargs...) = Object(Cuboid(), args...; kwargs...)


# special case: cubes

"""
    cube(x, y, z, w, v=1) (cube of width `w` centered at `(x,y,z)`)
    cube((x,y,z), w, v=1) ditto
    cube(w, v=1) centered at origin
Construct cubes as special cases of `Cuboid`.
"""
cube(cx::RealU, cy::RealU, cz::RealU, w::RealU, v::Number = 1) =
    cuboid(cx, cy, cz, w, w, w, 0, 0, v)
cube(center::NTuple{3,RealU}, w::RealU, v::Number = 1) =
    cube(center..., w, v)
cube(w::RealU, v::Number = 1) = cube((zero(w), zero(w), zero(w)), v)


# methods


volume1(::Cuboid) = 1 # volume of unit cube

ℓmax1(::Cuboid) = √3 # max line integral through unit cuboid

ℓmax(ob::Object3d{Cuboid}) = sqrt(sum(abs2, ob.width))


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
    u::Ru,
    v::Rv,
    ϕ::RealU, # azim
    θ::RealU; # polar
) where {Ru <: Real, Rv <: Real}
    T = promote_type(Ru, Rv, Float32)

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
    return T(ℓ)::T
end


# spectrum

"""
    spectrum1(::Object3d{Cuboid}, (kx,ky,kz))
Spectrum of unit cube at `(kx,ky,kz)`,
for unitless spatial frequency coordinates.
"""
function spectrum1(::Object3d{Cuboid}, kxyz::NTuple{3,Real})
    return prod(sinc, kxyz)
end
