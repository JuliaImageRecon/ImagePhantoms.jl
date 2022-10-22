#=
object3.jl
Utilities for 3D objects
=#

using LazyGrids: ndgrid
export phantom, radon, spectrum


volume(ob::Object3d{S}) where S = volume1(S()) * prod(ob.width)

ℓmax(ob::Object3d{S}) where S = ℓmax1(S()) * maximum(ob.width)


# rotate

"""
    rotate(ob::Object3d, (α,β))
    rotate(ob::Object3d, α, β=0)
Rotate a 3D object.
"""
rotate(ob::Object3d{S}, θ::NTuple{3,RealU}) where S =
    Object(S(), ob.center, ob.width, ob.angle .+ θ, ob.value)
rotate(ob::Object3d, α::RealU, β::RealU=0, γ::RealU=0) = rotate(ob, (α,β,γ))


"""
    coords(object::Object3d, x::RealU, y::RealU, y::RealU)
Put coordinates `(x,y,z)` in canonical axes associated with `object`.
"""
function coords(ob::Object3d, x::RealU, y::RealU, z::RealU)
    xyz = Rxyz_inv(x - ob.center[1], y - ob.center[2], z - ob.center[3],
        ob.sin..., ob.cos...)
    return xyz ./ ob.width # unitless
end


# phantom images


"""
    phantom(ob::Object3d)::Function
Return function of `(x,y,z)`
that user can sample at any locations
to make a 3D phantom image
for a single 3D object.
"""
function phantom(ob::Object3d)
    # apply coordinate transformation and value:
    return (x,y,z) -> ob.value * phantom1(ob, coords(ob, x, y, z))
end


"""
    phantom(oa::Array{<:Object3d})::Function
Return function of `(x,y,z)`
that user can sample at any locations
to make a 3D phantom image
for an `Array` of 3D objects.
"""
function phantom(oa::Array{<:Object3d})
    return (x,y,z) -> sum(ob -> phantom(ob)(x,y,z), oa)
end


"""
    phantom(x::Vector, y::Vector, z::Vector, oa::Array{<:Object3d})
Return a 3D digital image of the phantom sampled at grid of `(x,y,z)` locations.
Returned 3D image size is `length(x) × length(y) × length(z)`.
"""
function phantom(
    x::AbstractVector,
    y::AbstractVector,
    z::AbstractVector,
    oa::Array{<:Object3d},
)
    return phantom(oa).(ndgrid(x,y,z)...)
end


"""
    image = phantom(x, y, z, oa::Array{<:Object3d}, oversample::Int; T)
Return a 3D digital image of the phantom sampled
at grid of `(x,y,z)` locations,
with over-sampling factor `oversample` and element type `T`.
"""
function phantom(
    x::AbstractVector,
    y::AbstractVector,
    z::AbstractVector,
    oa::Array{<:Object3d{S,V}},
    oversample::Int;
    T::DataType = eltype(oneunit(V) * 1f0), # at least Float32
) where {S, V <: Number}

    oversample < 1 && throw(ArgumentError("oversample $oversample"))
    dx = x[2] - x[1]
    dy = y[2] - y[1]
    dz = z[2] - z[1]
    all(≈(dx), diff(x)) || throw("oversample requires uniform x")
    all(≈(dy), diff(y)) || throw("oversample requires uniform y")
    all(≈(dz), diff(z)) || throw("oversample requires uniform z")
    tmp = ((1:oversample) .- (oversample+1)/2) / oversample
    o3 = oversample^3
    ophantom = ob -> (x,y,z) ->
        T(sum(phantom(ob).(ndgrid(x .+ dx*tmp, y .+ dy*tmp, z .+ dz*tmp)...)) / o3)
    out = sum(ob -> ophantom(ob).(ndgrid(x,y,z)...), oa)
    return out::Array{T, 3}
end


# radon transform

# shift property of the X-ray transform; ϕ=azimuth θ=polar
function xray_shift(
    u::RealU, v::RealU, ϕ::RealU, θ::RealU, # projection coordinates
    cx::RealU, cy::RealU, cz::RealU, # object center
)
    (sϕ, cϕ) = sincos(ϕ)
    (sθ, cθ) = sincos(θ)
    ushift = cx * cϕ + cy * sϕ
    vshift = (cx * sϕ - cy * cϕ) * sθ + cz * cθ
    return (u - ushift, v - vshift, ϕ, θ)
end

# rotation property (only axial for now)
function xray_rotate(
    Δu::RealU, Δv::RealU, ϕ::RealU, θ::RealU, # projection coordinates
    Φazim::RealU, # object angles
    Θpolar::RealU,
    ψ::RealU,
)
    iszero(Θpolar) || throw("nonzero polar object angle not done")
    iszero(ψ) || throw("nonzero ψ object angle not done for radon") # todo
    return (Δu, Δv, ϕ - Φazim, θ)
end

# affine scaling property of the X-ray transform
# p,topo,prop,affine
function xray_scale(
    u::RealU, v::RealU, ϕ::RealU, θ::RealU, # projection coordinates
    wx::RealU, wy::RealU, wz::RealU, # object width
)
    (sϕ, cϕ) = sincos(ϕ)
    (sθ, cθ) = sincos(θ)
    e1old = [cϕ, sϕ, 0]
    e2old = [-sϕ*cθ, cϕ*cθ, sθ]
    e3old = [sϕ*sθ, -cϕ*sθ, cθ]
    width = [wx, wy, wz]
    r = 1 / sqrt( sum(abs2, e2old ./ width) ) # would be radius for sphere
    ϕp = atan(wy * sϕ, wx * cϕ) # ϕ'
    θp = asin(sθ * r / wz) # θ'
    tmp = (e1old * u + e3old * v) ./ width
    (sϕp, cϕp) = sincos(ϕp)
    (sθp, cθp) = sincos(θp)
    e1new = [cϕp, sϕp, 0]
    e3new = [sϕp*sθp, -cϕp*sθp, cθp]
    up = e1new' * tmp
    vp = e3new' * tmp
    return (r, up, vp, ϕp, θp)
end


# interface to xray1 after applying shift, rotate, scale properties
function _xray(
    type::AbstractShape{3},
    center::NTuple{3,RealU},
    width::NTuple{3,RealU},
    angle::NTuple{3,RealU},
    u::RealU, v::RealU, ϕ::RealU, θ::RealU,
)
    Δu, Δv, ϕ, θ = xray_shift(u, v, ϕ, θ, center...)
    ur, vr, ϕr, θr, = xray_rotate(Δu, Δv, ϕ, θ, angle...)
    scale, up, vp, ϕp, θp = xray_scale(ur, vr, ϕr, θr, width...)
    return scale * xray1(type, up, vp, ϕp, θp)
end


# this gateway seems to help type inference
function _radon(ob::Object3d{S}, u::RealU, v::RealU, ϕ::RealU, θ::RealU) where S
    T = radon_type(ob)
    return T(ob.value * _xray(S(), ob.center, ob.width, ob.angle, u, v, ϕ, θ))::T
end


"""
    radon(ob::Object3d)::Function
Return function of `(u,v,ϕ,θ)` that user can sample
at any projection coordinates
to make projection views of a 3D object.

The coordinate system used here is such that `ϕ=0` corresponds to
line integrals along the ``y`` axis
for an object ``f(x,y,z)``.
Then as `ϕ` increases, the line integrals rotate counter-clockwise.
"""
function radon(ob::Object3d{S}) where S
    return (u::RealU, v::RealU, ϕ::RealU, θ::RealU) -> _radon(ob, u, v, ϕ, θ)
end


"""
    radon(oa::Array{<:Object3d})::Function
Return function of `(u,v,ϕ,θ)` that user can sample
at any projection coordinates
to make projection views
for an array of 3D objects.
"""
function radon(oa::Array{<:Object3d})
    return (u::RealU, v::RealU, ϕ::RealU, θ::RealU) -> sum(ob -> radon(ob)(u,v,ϕ,θ), oa)
end


"""
    radon(u:Vector, v:Vector, ϕ:Vector, θ:Vector, oa::Array{<:Object3d})
Return parallel-beam projections
sampled at grid of `(u,v,ϕ,θ)` locations.
Returned array size is `length(u) × length(v) × length(ϕ) × length(θ)`.
"""
function radon(
    u::AbstractVector,
    v::AbstractVector,
    ϕ::AbstractVector,
    θ::AbstractVector,
    oa::Array{<:Object3d},
)
    return radon(oa).(ndgrid(u, v, ϕ, θ)...)
end


"""
    radon(u:Vector, v:Vector, ϕ:RealU, θ:RealU, oa::Array{<:Object3d})
Return parallel-beam projection view sampled at grid of `(u,v)` locations
for a given `(ϕ,θ)` pair.
Returned array size is `length(u) × length(v)`.
"""
function radon(
    u::AbstractVector,
    v::AbstractVector,
    ϕ::RealU,
    θ::RealU,
    oa::Array{<:Object3d},
)
    return radon(oa).(ndgrid(u, v)..., ϕ, θ)
end


# spectra


# apply rotate, translate, and scale properties of 3D Fourier transform
function _spectrum(ob::Object3d, fx, fy, fz, cx, cy, cz, rx, ry, rz) #, Φ, Θ, ψ)
    # rotate then translate:
#   (kx, ky, kz) = Rxyz_inv(fx, fy, fz, Φ, Θ, ψ)
    (kx, ky, kz) = Rxyz_inv(fx, fy, fz, ob.sin..., ob.cos...)
    return rx * ry * rz * cispi(-2*(fx*cx + fy*cy + fz*cz)) *
        spectrum1(ob, (rx*kx, ry*ky, rz*kz))
end


"""
    spectrum(ob::Object3d)::Function
Return function of `(fx,fy,fz)` that user can sample
at any spatial frequency locations
to evaluate the spectrum (3D Fourier transform)
of a single 3D object.
Units of spatial frequencies should be the reciprocal
of the units defining the object.
"""
function spectrum(ob::Object3d)
    return (fx,fy,fz) -> ob.value *
        _spectrum(ob, fx, fy, fz, ob.center..., ob.width...) #, ob.angle...)
end


"""
    spectrum(oa::Array{<:Object3d})::Function
Return function `kspace(fx,fy,fz)` that user can sample
at any spatial frequency locations
to evaluate the spectrum (3D Fourier transform)
for an `Array` of 3D objects.
"""
function spectrum(oa::Array{<:Object3d})
    return (fx,fy,fz) -> sum(ob -> spectrum(ob)(fx,fy,fz), oa)
end


"""
    spectrum(fx::Vector, fy::Vector, fz::Vector, oa::Array{<:Object3d})
Return 3D k-space array
sampled at grid of `(fx,fy,fz)` locations.
Returned 3D array size is `length(fx) × length(fy) × length(fz)`.
"""
function spectrum(
    fx::AbstractVector,
    fy::AbstractVector,
    fz::AbstractVector,
    oa::Array{<:Object3d},
)
    return spectrum(oa).(ndgrid(fx,fy,fz)...)
end
