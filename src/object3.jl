#=
object3.jl
Utilities for 3D objects, cf object.jl
=#

using LazyGrids: ndgrid

# methods for phantoms: an array of objects

"""
    image = phantom(oa::Array{<:Object3d})::Function
Return function `image(x,y,z)` that user can sample at any `(x,y,z)` locations
to make a phantom image.
"""
function phantom(oa::Array{<:Object3d})
    return (x,y,z) -> sum(ob -> phantom(ob)(x,y,z), oa)
end

"""
    image = phantom(x, y, oa::Array{<:Object3d}, oversample::Int; T)
Return a digital image of the phantom sampled at `(x,y,z)` locations,
with over-sampling factor `oversample` and element type `T`.
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
    ophantom = ob ->
       (x,y,z) -> T(sum(phantom(ob).(ndgrid(x .+ dx*tmp, y .+ dy*tmp, z .+ dz*tmp)...)) / abs2(oversample))
    return sum(ob -> ophantom(ob).(ndgrid(x,y,z)...), oa)
end

"""
    image = phantom(x, y, z, oa::Array{<:Object3d})
Return a digital image of the phantom sampled at `(x,y,z)` locations.
"""
function phantom(
    x::AbstractVector,
    y::AbstractVector,
    z::AbstractVector,
    oa::Array{<:Object3d},
)

#   return sum(ob -> phantom(ob).(ndgrid(x,y,z)...), oa) # todo
    return phantom(ndgrid(x,y,z)..., oa)
end

function phantom(
    xx::AbstractArray,
    yy::AbstractArray,
    zz::AbstractArray,
    oa::Array{<:Object3d},
)
    return sum(ob -> phantom(ob).(xx,yy,zz), oa)
end


"""
    sino = radon(oa::Array{<:Object3d})::Function
Return function `sino(u,v,ϕ,θ)` that user can sample at any `(u,v,ϕ,θ)` locations
to make projection views of a 3D object.

The coordinate system used here is such that `ϕ=0` corresponds to
line integrals along the ``y`` axis for an object ``f(x,y,z)``.
Then as `ϕ` increases, the line integrals rotate counter-clockwise.
"""
function radon(oa::Array{<:Object3d})
    return (u,v,ϕ,θ) -> sum(ob -> radon(ob)(u,v,ϕ,θ), oa)
end

"""
    sino = radon(u, v, ϕ, θ, oa::Array{<:Object3d})
Return parallel-beam projections `sino` sampled at given `(u,v,ϕ,θ)` locations.
"""
function radon(
    u::AbstractVector,
    v::AbstractVector,
    ϕ::AbstractVector,
    θ::AbstractVector,
    oa::Array{<:Object3d},
)
#   return sum(ob -> radon(ob).(ndgrid(u, v, ϕ, θ)...), oa) # todo
    return radon(ndgrid(u, v, ϕ, θ)..., oa)
end

function radon(
    u::AbstractArray,
    v::AbstractArray,
    ϕ::AbstractArray,
    θ::AbstractArray,
    oa::Array{<:Object3d},
)
    return sum(ob -> radon(ob).(u,v,ϕ,θ), oa)
end


"""
    kspace = spectrum(oa::Array{<:Object3d})::Function
Return function `kspace(fx,fy,fz)` that user can sample at any `(fx,fy,fz)` locations
to make phantom 3D k-space data.
"""
function spectrum(oa::Array{<:Object3d})
    return (fx,fy,fz) -> sum(ob -> spectrum(ob)(fx,fy,fz), oa)
end

"""
    kspace = spectrum(fx, fy, oa::Array{<:Object3d})::Function
Return k-space array `kspace` sampled at given `(fx,fy,fz)` locations.
"""
function spectrum(
    fx::AbstractVector,
    fy::AbstractVector,
    fz::AbstractVector,
    oa::Array{<:Object3d},
)
#   return sum(ob -> spectrum(ob).(fx,fy'), oa) # todo
    return spectrum(ndgrid(fx,fy,fz)..., oa)
end

function spectrum(
    fx::AbstractArray,
    fy::AbstractArray,
    fz::AbstractArray,
    oa::Array{<:Object3d},
)
    return sum(ob -> spectrum(ob).(fx,fy,fz), oa)
end


# helpers


function rotate3d(x::RealU, y::RealU, z::RealU, ϕ::RealU, θ::RealU)
    θ == 0 || throw(ArgumentError("θ ≂̸ 0 unsupported currently"))
    (s, c) = sincos(ϕ)
    return (c * x + s * y, -s * x + c * y, z)
end

rotate3d(xyz::NTuple{3,RealU}, ϕ::RealU, θ::RealU) = rotate3d(xyz..., ϕ, θ)
