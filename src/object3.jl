#=
object3.jl
Utilities for 3D objects, cf object.jl
=#

using LazyGrids: ndgrid


# phantom images


"""
    phantom(ob::Object3d)::Function
Return function of `(x,y,z)` that user can sample at any `(x,y,z)` locations
to make a 3D phantom image for a single 3D object.
"""
function phantom(ob::Object3d)
#   apply coordinate transformation and value:
    return (x,y,z) -> ob.value * phantom1(ob, coords(ob, x, y, z))
end


"""
    phantom(oa::Array{<:Object3d})::Function
Return function of `(x,y,z)` that user can sample at any `(x,y,z)` locations
to make a 3D phantom image for an `Array` of 3D objects.
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
Return a 3D digital image of the phantom sampled at grid of `(x,y,z)` locations,
with over-sampling factor `oversample` and returned element type `T`.
"""
function phantom(
    x::AbstractVector,
    y::AbstractVector,
    z::AbstractVector,
    oa::Array{<:Object3d},
    oversample::Int;
    T::DataType = promote_type(eltype.(oa)..., Float32),
)

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
    return sum(ob -> ophantom(ob).(ndgrid(x,y,z)...), oa)
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

# affine scaling property of the X-ray transform; ϕ=azimuth θ=polar
# p,topo,prop,affine
function xray_scale(
    u::RealU, v::RealU, ϕ::RealU, θ::RealU, # projection coordinates
    wx::RealU, wy::RealU, wz::RealU, # object width
)
    (sϕ, cϕ) = sincos(ϕ)
    ϕp = atan(wy * sϕ, wx * cϕ) # ϕ'
    (sθ, cθ) = sincos(θ)
    denom = sqrt((wx * cϕ)^2 + (wy * sϕ)^2)
    θp = atan(u * v * sθ, v * denom) # θ'
    vp = (wy/wx - wx/wy) * sϕ * cϕ * u + cos(θp)/cos(θ) * v / c
    return (u / denom, vp, ϕp, θp)
end


function _radon(
    type::AbstractShape3,
    u::Real, v::Real, ϕ::RealU, θ::RealU,
    xang::RealU, zang::RealU,
)
throw("todo") # need a rotation property too
    return xray1(type, u, v, ϕ, θ)
end


# todo: not finished
"""
    radon(ob::Object3d)::Function
Return function of `(u,v,ϕ,θ)` that user can sample
at any projection coordinates
to make projection views of a 3D object.

The coordinate system used here is such that `ϕ=0` corresponds to
line integrals along the ``y`` axis for an object ``f(x,y,z)``.
Then as `ϕ` increases, the line integrals rotate counter-clockwise.
"""
function radon(ob::Object3d)
    return (u,v,ϕ,θ) -> ob.value *
        _radon(
         xray_scale(
          xray_shift(u, v, ϕ, θ, ob.center...)...,
          ob.width...,
         ),
         ob.angles...,
        )
end


"""
    radon(oa::Array{<:Object3d})::Function
Return function of `(u,v,ϕ,θ)` that user can sample
at any projection coordinates
to make projection views of an array of 3D objects.
"""
function radon(oa::Array{<:Object3d})
    return (u,v,ϕ,θ) -> sum(ob -> radon(ob)(u,v,ϕ,θ), oa)
end


"""
    radon(u:Vector, v:Vector, ϕ:Vector, θ:Vector, oa::Array{<:Object3d})
Return parallel-beam projections sampled at grid of `(u,v,ϕ,θ)` locations.
Returned array size is `length(u) × length(v) × length(ϕ) × length(θ)`.
"""
function radon(
    u::AbstractVector,
    v::AbstractVector,
    ϕ::AbstractVector,
    θ::AbstractVector,
    oa::Array{<:Object3d},
)
#   return sum(ob -> radon(ob).(ndgrid(u, v, ϕ, θ)...), oa) # todo
    return radon(oa).(ndgrid(u, v, ϕ, θ)...)
end

#=
function radon(
    u::AbstractArray,
    v::AbstractArray,
    ϕ::AbstractArray,
    θ::AbstractArray,
    oa::Array{<:Object3d},
)
    return sum(ob -> radon(ob).(u,v,ϕ,θ), oa)
end
=#


# spectra


# apply rotate, translate, and scale properties of 3D Fourier transform
function _spectrum(ob::Object3d, fx, fy, fz, cx, cy, cz, rx, ry, rz, Φ, Θ)
    (kx, ky, kz) = rotate3d(fx, fy, fz, Φ, Θ) # rotate then translate
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
        _spectrum(ob, fx, fy, fz, ob.center..., ob.width..., ob.angle...)
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

#=
function spectrum(
    fx::AbstractArray,
    fy::AbstractArray,
    fz::AbstractArray,
    oa::Array{<:Object3d},
)
    return sum(ob -> spectrum(ob).(fx,fy,fz), oa)
end
=#


# helpers


function rotate3d(x::RealU, y::RealU, z::RealU, ϕ::RealU, θ::RealU)
    θ == 0 || throw(ArgumentError("θ ≂̸ 0 unsupported currently"))
    (s, c) = sincos(ϕ)
    return (c * x + s * y, -s * x + c * y, z)
end

rotate3d(xyz::NTuple{3,RealU}, ϕ::RealU, θ::RealU) = rotate3d(xyz..., ϕ, θ)
