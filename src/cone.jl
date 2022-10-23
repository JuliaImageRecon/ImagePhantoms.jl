#=
cone.jl
Right cone
Base object is a cone with unit radius and unit height,
with base at origin pointing up along z.
=#


export Cone, cone


"""
    Cone <: AbstractShape{3}
"""
struct Cone <: AbstractShape{3} end


# constructor


"""
    cone(cx, cy, cz, wx, wy, wz, Φ, Θ, value::Number)
    cone(center::NTuple{3,RealU}, width::NTuple{3,RealU}, angle::NTuple{3,RealU}, v)
Construct `Object{Cone}` from parameters;
here `width` is the base radii for `wx` and `wy` and `wz` is the height.
"""
cone(args... ; kwargs...) = Object(Cone(), args...; kwargs...)


# methods


volume1(::Cone) = π/3 # volume of unit cone

ℓmax1(::Cone) = 2 # max line integral through unit cone

function ℓmax(ob::Object3d{Cone})
    rmax = maximum(ob.width[1:2])
    return max(2 * rmax, sqrt(rmax^2 + ob.width[3]^2))
end


"""
    phantom1(ob::Object3d{Cone}, (x,y,z))
Evaluate unit cone at `(x,y,z)`,
for unitless coordinates.
"""
function phantom1(ob::Object3d{Cone}, xyz::NTuple{3,Real})
    z = xyz[3]
    r = sqrt(sum(abs2, xyz[1:2]))
    return (0 ≤ z ≤ 1) && r ≤ 1 - z
end


# radon


# x-ray transform (line integral) of unit cone
# `u,v` should be unitless
#=
function xray1(
    ::Cone,
    u::Ru,
    v::Rv,
    ϕ::RealU, # azim
    θ::RealU; # polar
) where {Ru <: Real, Rv <: Real}
    T = promote_type(Ru, Rv, Float32)

#   (sϕ, cϕ) = sincos(ϕ)
#   (sθ, cθ) = sincos(θ)

    throw("todo: cone radon not done")
end
=#


# spectrum

#=
"""
    spectrum1(::Object3d{Cone}, (kx,ky,kz))
Spectrum of unit cone at `(kx,ky,kz)`,
for unitless spatial frequency coordinates.
"""
function spectrum1(::Object3d{Cone}, kxyz::NTuple{3,Real})
    throw("todo: cone spectrum not done")
end
=#
