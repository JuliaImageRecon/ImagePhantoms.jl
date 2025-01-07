#=
object.jl
Types for making (often) piece-wise constant phantoms
out of elementary objects having shapes like
ellipses, rectangles, and gaussian bumps.
=#

export AbstractObject
export Object, Object2d, Object3d
export rotate, scale, translate
export phantom, radon, spectrum

atuple(x::Any, n::Int) = ntuple(i -> x, n)

abstract type AbstractObject end


"""
    Object{S, D, V, ...}(center, width, angle, value, param) <: AbstractObject
General container for 2D and 3D objects for defining image phantoms.

* `center::NTuple{D,C}` coordinates of "center" of this object
* `width::NTuple{W,C}` "width" along axis; usually `W=D` or `W=1`
  (e.g., FWHM for Gauss, radii for Ellipse, radius of Circle)
* `angle::NTuple{D-1,A}` angle of x' axis relative to x axis,
  in radians (or with units)
* `value::V` "intensity" value for this object
* `param` optional additional shape parameters (often `nothing`)

Example

`jldoctest
julia> Object(Ellipse(), (0,0), (1,2), 0.0, 1//2, nothing)

Object2d{Ellipse, Rational{Int64}, 2, Int64, Float64, Nothing} (S, D, V, ...)
 center::Tuple{Int64, Int64} (0, 0)
 width::Tuple{Int64, Int64} (1, 2)
 angle::Tuple{Float64} (0.0,)
 value::Rational{Int64} 1//2
 param::Nothing nothing
`
"""
struct Object{S, D, V, W, C, A, P} <: AbstractObject
    shape::S
    "x,y center coordinates"
    center::NTuple{D,C}
    "'width' along x',y' axes (FWHM for Gauss, radii for Ellipse)"
    width::NTuple{W,C}
    "angle of x' axis relative to x axis, in radians (or with units)"
    angle::NTuple{Da,A} where Da
    "'intensity' value for this shape"
    value::V
    "optional additional shape parameters"
    param::P # often `nothing`

    """
        shape = Object(shape, center, width, angle, value, param)
    Top-level shape constructor where `angle` must be tuple.
    """
    function Object{S}(
        shape::S,
        center::NTuple{D,C},
        width::NTuple{W,C},
        angle::NTuple{Da,A},
        value::V,
        param::P,
    ) where {S <: AbstractShape, D, V <: Number, W, C <: RealU, Da, A <: RealU, P}
        1 ≤ Da == D-1 || throw(ArgumentError("Dϕ=$Dϕ != D-1, where D=$D"))
        all(width .> zero(eltype(width))) || throw(ArgumentError("widths must be positive"))
        new{S,D,V,W,C,A,P}(shape, center, width, angle, value, param)
    end
end


"""
    Object2d = Object{S,2} where S <: AbstractObject
For 2D objects
"""
const Object2d = Object{S,2} where S

"""
    Object3d = Object{S,3} where S <: AbstractObject
For 3D objects
"""
const Object3d = Object{S,3} where S


# Constructors

"""
    ob = Object(shape, center, width, angle, value, param)
General outer object constructor where `angle` is a tuple, including 3D case.
"""
function Object(
    shape::S,
    center::NTuple{D,C},
    width::NTuple{W,C} where W = atuple(zero(eltype(center)), length(center)),
    angle::NTuple{Da,RealU} where Da = atuple(0, length(center)-1),
    value::Number = 1,
    param::Any = nothing,
) where {S <: AbstractShape, D, C <: RealU}
    Object{S}(shape, center, width, angle, value, param)
end


"""
    ob = Object(shape, center, width, angle=0, value=1, param=nothing)
2D object constructor.
"""
function Object(
    shape::AbstractShape2,
    center::NTuple{2,C},
    width::NTuple{W,C} where W,
    angle::RealU = 0,
    value::Number = 1,
    param::Any = nothing,
) where {C <: RealU}
    Object(shape, center, width, (angle,), value, param)
end

"""
    ob = Object(shape ; center=(0,0), width=(1,1), angle=0, value=1, param=nothing)
2D object constructor using keywords.
"""
function Object(
    shape::AbstractShape2 ;
    center::NTuple{2,C} = (0,0),
    width::NTuple{W,C} where W = atuple(one(eltype(center)), 2),
    angle::RealU = 0,
    value::Number = 1,
    param = nothing,
) where {C <: RealU}
    Object(shape, center, width, angle, value, param)
end

#=
# not needed due to general case above
"""
    ob = Object(shape, center, width, angle=(0,0), value=1, param=nothing)
3D object constructor.
"""
function Object(
    shape::Shape3d,
    center::NTuple{3,C},
    width::NTuple{W,C} where W,
    angle::NTuple{2,RealU} = (0,0),
    value::Number = 1,
    param::Any = nothing,
) where {C <: RealU}
    Object(shape, center, width, angle, value, param)
end
=#

"""
    ob = Object(shape ; center=(0,0,0), width=(1,1,1), angle=(0,0), value=1, param=nothing)
3D object constructor using keywords.
"""
function Object(
    shape::AbstractShape3 ;
    center::NTuple{3,C} = (0,0,0),
    width::NTuple{W,C} = atuple(one(eltype(center)), 3),
    angle::NTuple{2,RealU} = (0,0),
    value::Number = 1,
    param = nothing,
) where {W, C <: RealU}
    Object(shape, center, width, angle, value, param)
end


#=
"""
    Object(shape ; x, y, rx, ry, ϕ=0, value=1 ; param=nothing)
2D object constructor without tuples.
"""
function Object(
    shape::AbstractShape,
    cx::C, cy::C,
    rx::C, ry::C,
    ϕ::RealU,
    value::Number = 1 ;
    param = nothing,
) where {C <: RealU}
# todo: might want to use promote
    Object(shape, (cx, cy), (rx, ry), ϕ, value, param)
end
=#


# Methods for objects

Base.eltype(::Object{S,D,V}) where {S,D,V} = V


"""
    show(io::IO, ::MIME"text/plain", ob::Object)
"""
function Base.show(io::IO, ::MIME"text/plain", ob::Object{S,D}) where {S,D}
    println(io, typeof(ob), " (S, D, V, ...)")
    for f in (:shape, :center, :width, :angle, :value, :param)
        p = getproperty(ob, f)
        println(io, " ", f, "::", typeof(p), " ", p)
    end
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
    scale(ob, atuple(factor,W))


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
    image = phantom(x, y, oa::Array{<:Object2d})
Return a digital image of the phantom sampled at `(x,y)` locations.
"""
function phantom(x::AbstractVector, y::AbstractVector, oa::Array{<:Object2d})
    return sum(ob -> phantom(ob).(x,y'), oa)
end

function phantom(xx::AbstractArray, yy::AbstractArray, oa::Array{<:Object2d})
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
function radon(r::AbstractVector, ϕ::AbstractVector, oa::Array{<:Object2d})
    return sum(ob -> radon(ob).(r,ϕ'), oa)
end

function radon(rr::AbstractArray, ϕϕ::AbstractArray, oa::Array{<:Object2d})
    return sum(ob -> radon(ob).(rr,ϕϕ), oa)
end


"""
    kspace = spectrum(oa::Array{<:Object2d})::Function
Return function `kspace(fx,fy)` that user can sample at any `(fx,fy)` locations
to make phantom 2D k-space data.
"""
function spectrum(oa::Array{<:Object2d})
    return (fx,fy) -> sum(ob -> spectrum(ob)(fx,fy), oa)
end

"""
    kspace = spectrum(fx, fy, oa::Array{<:Object2d})::Function
Return k-space matrix `kspace` sampled at given `(fx,fy)` locations.
"""
function spectrum(fx::AbstractVector, fy::AbstractVector, oa::Array{<:Object2d})
    return sum(ob -> spectrum(ob).(fx,fy'), oa)
end

function spectrum(fx::AbstractArray, fy::AbstractArray, oa::Array{<:Object2d})
    return sum(ob -> spectrum(ob).(fx,fy), oa)
end


# helpers


function rotate2d(x::RealU, y::RealU, θ::RealU)
    (s, c) = sincos(θ)
    return (c * x + s * y, -s * x + c * y)
end

rotate2d(xy::NTuple{2,RealU}, θ::RealU) = rotate2d(xy..., θ)
