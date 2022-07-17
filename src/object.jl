#=
object.jl
Types for making (often) piece-wise constant phantoms
out of elementary objects having shapes like
ellipses, rectangles, and gaussian bumps.
=#

export AbstractObject
export Object, Object2d, Object3d
#export rotate, scale, translate # no: to avoid conflicts with Plots
export phantom, radon, spectrum

_tuple(x::Any, n::Int) = ntuple(i -> x, n)

abstract type AbstractObject end


"""
    Object{S, D, V, ...}(center, width, angle, value, param) <: AbstractObject
General container for 2D and 3D objects for defining image phantoms.

* `center::NTuple{D,C}` coordinates of "center" of this object
* `width::NTuple{D,C}` "width" along each axis
  (e.g., FWHM for Gauss, radii for Ellipse)
* `angle::NTuple{D-1,A}` angle of x' axis relative to x axis,
  in radians (or with units)
* `value::V` "intensity" value for this object
* `param` optional additional shape parameters (often `nothing`)

# Example

```jldoctest
julia> Object(Ellipse(), (0,0), (1,2), 0.0, 1//2, nothing)
Object2d{Ellipse, Rational{Int64}, 2, Int64, Float64, Nothing} (S, D, V, ...)
 shape::Ellipse Ellipse()
 center::NTuple{2,Int64} (0, 0)
 width::NTuple{2,Int64} (1, 2)
 angle::Tuple{Float64} (0.0,)
 value::Rational{Int64} 1//2
 param::Nothing nothing
```
"""
struct Object{S, D, V, C, A, P} <: AbstractObject
    shape::S
    "x,y center coordinates"
    center::NTuple{D,C}
    "'width' along x',y' axes (FWHM for Gauss, radii for Ellipse)"
    width::NTuple{D,C}
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
        width::NTuple{D,C},
        angle::NTuple{Da,A},
        value::V,
        param::P,
    ) where {S <: AbstractShape, D, V <: Number, C <: RealU, Da, A <: RealU, P}
        1 ≤ Da == D-1 || throw(ArgumentError("Da=$Da != D-1, where D=$D"))
        all(width .> zero(eltype(width))) || throw(ArgumentError("widths must be positive"))
        new{S,D,V,C,A,P}(shape, center, width, angle, value, param)
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
    width::NTuple{D,C} = _tuple(zero(eltype(center)), length(center)),
    angle::NTuple{Da,RealU} where Da = _tuple(0, length(center)-1),
    value::Number = 1,
    param::Any = nothing,
) where {S <: AbstractShape, D, C <: RealU}
    Object{S}(shape, center, width, promote(angle...), value, param)
end


"""
    ob = Object(shape, center, width, angle=0, value=1, param=nothing)
2D object constructor using tuples.
"""
function Object(
    shape::AbstractShape{2},
    center::NTuple{2,C},
    width::NTuple{2,C},
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
    shape::AbstractShape{2} ;
    center::NTuple{2,C} = (0,0),
    width::NTuple{2,C} = _tuple(one(eltype(center)), 2),
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
3D object constructor using tuples.
"""
function Object(
    shape::Shape{3},
    center::NTuple{3,C},
    width::NTuple{3,C},
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
    shape::AbstractShape{3} ;
    center::NTuple{3,C} = (0,0,0),
    width::NTuple{3,C} = _tuple(one(eltype(center)), 3),
    angle::NTuple{2,RealU} = (0,0),
    value::Number = 1,
    param = nothing,
) where {C <: RealU}
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


"""
    Object(ob::Object ; center, width, angle, value, param)
Make a copy of Object `ob`, optionally modifying some values.
"""
function Object(ob::Object{S,D} ;
    center::NTuple{D} = ob.center,
    width::NTuple{D} = ob.width,
    angle::NTuple{Da,<:RealU} = ob.angle,
    value::Number = ob.value,
    param = ob.param,
) where {S, D, Da}
    Da == D-1 || throw(ArgumentError("Da=$Da != D-1, where D=$D"))
    Object(ob.shape, center, width, angle, value, param)
end


# Methods for objects

Base.eltype(::Object{S,D,V}) where {S,D,V} = V
Base.ndims(::Object{S,D}) where {S,D} = D


"""
    show(io::IO, ::MIME"text/plain", ob::Object)
"""
function Base.show(io::IO, ::MIME"text/plain", ob::Object{S,D}) where {S,D}
    println(io, typeof(ob), " (S, D, V, ...)")
    for f in (:shape, :center, :width, :angle, :value, :param)
        p = getproperty(ob, f)
        t = typeof(p)
        t = t == NTuple{D,eltype(t)} ? "NTuple{$D,$(eltype(t))}" : "$t"
        println(io, " ", f, "::", t, " ", p)
    end
end


# scale width (tricky because length(width) < D is possible
"""
    scale(ob::Object, factor::RealU)
    scale(ob::Object, factor::NTuple{D,RealU})
Scale the width(s) by `factor`.
"""
scale(ob::Object{S,D}, factor::NTuple{D,RealU}) where {S,D} =
    Object(ob.shape, ob.center, ob.width .* factor, ob.angle, ob.value, ob.param)
scale(ob::Object{S,D}, factor::RealU) where {S,D} =
    scale(ob, _tuple(factor,D))


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
