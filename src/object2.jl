#=
object2.jl
Utilities for 2D objects
=#

#export rotate, scale, translate # no: to avoid conflicts with Plots
export phantom, radon, spectrum

_tuple(x::Any, n::Int) = ntuple(i -> x, n)



end


# rotate
"""
    rotate(ob::Object2d, θ::RealU)
In-plane rotation of a 2D object.
"""
rotate(ob::Object2d, θ::RealU) =
    Object(ob.shape, ob.center, ob.width, ob.angle .+ (θ,), ob.value, ob.param)
#rotate(ob::Object3d, θ::RealU) =
#Object(ob.shape, ob.center, ob.width, ob.angle .+ (θ,0), ob.value, ob.param)


"""
    rotate(ob::Object3d, (α,β))
    rotate(ob::Object3d, α, β=0)
Rotation of a 3D object.
"""
rotate(ob::Object3d, θ::NTuple{2,RealU}) =
    Object(ob.shape, ob.center, ob.width, ob.angle .+ θ, ob.value, ob.param)
rotate(ob::Object3d, α::RealU, β::RealU=0) = rotate(ob, (α,β))


# scale width (tricky because length(width) < D is possible
"""
    scale(ob::Object, factor::RealU)
    scale(ob::Object, factor::NTuple{W,RealU})
Scale the width(s) by `factor`.
"""
scale(ob::Object{S,D,V,W}, factor::NTuple{W,RealU}) where {S,D,V,W} =
    Object(ob.shape, ob.center, ob.width .* factor, ob.angle, ob.value, ob.param)
scale(ob::Object{S,D,V,W}, factor::RealU) where {S,D,V,W} =
    scale(ob, _tuple(factor,W))


# scale value
"""
    (*)(ob::Object, x::Number)
    (*)(x::Number, ob::Object)
Scale object `value` by `x`.
"""
Base.:(*)(ob::AbstractObject, x::Number) =
    Object(ob.shape, ob.center, ob.width, ob.angle, ob.value * x, ob.param)
Base.:(/)(ob::Object, x::Number) = ob * (1 / x)
Base.:(*)(x::Number, ob::Object) = ob * x


# translate
"""
    translate(ob::Object, shift::NTuple{D,RealU})
    translate(ob::Object2d, xshift, yshift)
    translate(ob::Object3d, xshift, yshift, zshift)
Translate the center coordinates of an object by `shift`
"""
translate(ob::Object{S,D}, shift::NTuple{D,RealU}) where {S, D} =
    Object(ob.shape, ob.center .+ shift, ob.width, ob.angle, ob.value, ob.param)
translate(ob::Object3d, x::RealU, y::RealU, z::RealU) = translate(ob, (x,y,z))
translate(ob::Object2d, x::RealU, y::RealU) = translate(ob, (x,y))


"""
    coords(object::Object2d, x::RealU, y::RealU)
Put coordinates `(x,y)` in canonical axes associated with `object`.
"""
function coords(ob::Object2d, x::RealU, y::RealU)
    (x, y) = rotate2d(x - ob.center[1], y - ob.center[2], ob.angle[1])
    return (x / ob.width[1], y / ob.width[2]) # unitless
end


"""
    coords(object::Object3d, x::RealU, y::RealU, y::RealU)
Put coordinates `(x,y,z)` in canonical axes associated with `object`.
"""
function coords(ob::Object3d, x::RealU, y::RealU, z::RealU)
    (x, y, z) = rotate3d(x - ob.center[1], y - ob.center[2], z - ob.center[3],
        ob.angle[1], ob.angle[2])
    return (x / ob.width[1], y / ob.width[2], z / ob.width[3]) # unitless
end


# methods for phantoms: an array of objects

"""
    image = phantom(oa::Array{<:Object2d})::Function
Return function `image(x,y)` that user can sample at any `(x,y)` locations
to make a phantom image.
"""
function phantom(oa::Array{<:Object2d})
    return (x,y) -> sum(ob -> phantom(ob)(x,y), oa)
end

"""
    image = phantom(x, y, oa::Array{<:Object2d}, oversample::Int; T)
Return a digital image of the phantom sampled at `(x,y)` locations,
with over-sampling factor `oversample` and element type `T`.
"""
function phantom(
    x::AbstractVector,
    y::AbstractVector,
    oa::Array{<:Object2d},
    oversample::Int;
    T::DataType = promote_type(eltype.(oa)..., Float32),
)
    oversample < 1 && throw(ArgumentError("oversample $oversample"))
    dx = x[2] - x[1]
    dy = y[2] - y[1]
    all(≈(dx), diff(x)) || throw("oversample requires uniform x")
    all(≈(dy), diff(y)) || throw("oversample requires uniform y")
    tmp = ((1:oversample) .- (oversample+1)/2) / oversample
    ophantom = ob ->
       (x,y) -> T(sum(phantom(ob).(x .+ dx*tmp, y .+ dy*tmp')) / abs2(oversample))
    return sum(ob -> ophantom(ob).(x,y'), oa)
end

"""
    image = phantom(x, y, oa::Array{<:Object2d})
Return a digital image of the phantom sampled at `(x,y)` locations.
"""
function phantom(
    x::AbstractVector,
    y::AbstractVector,
    oa::Array{<:Object2d},
)

    return sum(ob -> phantom(ob).(x,y'), oa)
end

function phantom(
    xx::AbstractArray,
    yy::AbstractArray,
    oa::Array{<:Object2d},
)
    return sum(ob -> phantom(ob).(xx,yy), oa)
end


"""
    sino = radon(oa::Array{<:Object2d})::Function
Return function `sino(r,ϕ)` that user can sample at any `(r,ϕ)` locations
to make a phantom 2D sinogram.

The coordinate system used here is such that `ϕ=0` corresponds to
line integrals along the ``y`` axis for an object ``f(x,y)``.
Then as `ϕ` increases, the line integrals rotate counter-clockwise.
"""
function radon(oa::Array{<:Object2d})
    return (r,ϕ) -> sum(ob -> radon(ob)(r,ϕ), oa)
end


"""
    sino = radon(r, ϕ, oa::Array{<:Object2d})
Return parallel-beam 2D sinogram `sino` sampled at given `(r,ϕ)` locations.
"""
function radon(
    r::AbstractVector,
    ϕ::AbstractVector,
    oa::Array{<:Object2d},
)
    return sum(ob -> radon(ob).(r,ϕ'), oa)
end


# spectra


# apply rotate, translate, and scale properties of 2D Fourier transform
function _spectrum(ob::Object2d, fx, fy, cx, cy, rx, ry, Φ)
    (kx, ky) = rotateld(fx, fy, Φ) # rotate then translate
    return rx * ry * cispi(-2*(fx*cx + fy*cy)) *
        spectrum1(ob, (rx*kx, ry*ky))
end


"""
    spectrum(ob::Object2d)::Function
Return function of `(fx,fy)` that user can sample
at any spatial frequency locations
to evaluate the spectrum (2D Fourier transform)
of a single 2D object.
Units of spatial frequencies should be the reciprocal
of the units defining the object.
"""
function spectrum(ob::Object2d)
    return (fx,fy) -> ob.value *
        _spectrum(ob, fx, fy, ob.center..., ob.width..., ob.angle)
end


"""
    spectrum(oa::Array{<:Object2d})::Function
Return function `kspace(fx,fy)` that user can sample
at any spatial frequency locations
to evaluate the spectrum (2D Fourier transform)
for an `Array` of 2D objects.
"""
function spectrum(oa::Array{<:Object2d})
    return (fx,fy) -> sum(ob -> spectrum(ob)(fx,fy), oa)
end


"""
    spectrum(fx:Vector, fy:Vector, oa::Array{<:Object2d})
Return 2D k-space array
sampled at grid of `(fx,fy)` locations.
Returned 2D array size is `length(fx) × length(fy)`.
"""
function spectrum(
    fx::AbstractVector,
    fy::AbstractVector,
    oa::Array{<:Object2d},
)
    return spectrum(oa).(fx,fy')
end


# helpers


function rotate2d(x::RealU, y::RealU, θ::RealU)
    (s, c) = sincos(θ)
    return (c * x + s * y, -s * x + c * y)
end

rotate2d(xy::NTuple{2,RealU}, θ::RealU) = rotate2d(xy..., θ)
