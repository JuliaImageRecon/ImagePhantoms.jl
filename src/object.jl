#=
object.jl
Types for making (often piece-wise constant) phantoms
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
    Object{S, D, V, ...}(center, width, angle, value) <: AbstractObject
General container for 2D and 3D objects for defining image phantoms.

* `center::NTuple{D,C}` coordinates of "center" of this object
* `width::NTuple{D,C}` "width" along each axis
  (e.g., FWHM for Gauss, radii for Ellipse)
* `angle::NTuple{D-1,A}` angle of x' axis relative to x axis,
  in radians (or with units)
* `value::V` "intensity" value for this object

# Example

```jldoctest
julia> Object(Ellipse(), (0,0), (1,2), 0.0, 1//2, nothing)
Object2d{Ellipse, Rational{Int64}, Int64, Float64, Nothing, 1} (S, D, V, ...)
 center::NTuple{2,Int64} (0, 0)
 width::NTuple{2,Int64} (1, 2)
 angle::Tuple{Float64} (0.0,)
 value::Rational{Int64} 1//2
```
"""
struct Object{S, D, V, C, A, Da} <: AbstractObject
    "x,y center coordinates"
    center::NTuple{D,C}
    "'width' along x',y' axes (FWHM for Gauss, radii for Ellipse)"
    width::NTuple{D,C}
    "angle of x' axis relative to x axis, in radians (or with units)"
    angle::NTuple{Da,A}
    "'intensity' value for this shape"
    value::V

    """
        Object{S}(center, width, angle, value)
    Inner constructor for `S <: AbstractShape`.
    The `center` and `width` tuples should have the same units
    (or should both be unitless).
    """
    function Object{S}(
        center::NTuple{D,RealU},
        width::NTuple{D,RealU},
        angle::Union{RealU, NTuple{Da,RealU} where Da},
        value::V,
    ) where {S <: AbstractShape, D, V <: Number}
        D == ndims(S()) || throw(ArgumentError("D=$D vs ndims(S)=$(ndims(S)) for S=$S"))
        if D == 2 && angle isa RealU
            angle = (angle,)
        end
        angle = promote(angle...)
        Da = length(angle)
        1 ≤ Da == D-1 || throw(ArgumentError("Da=$Da != D-1, where D=$D"))

        all(width .> zero(eltype(width))) || throw(ArgumentError("widths must be positive"))

        C = promote_type(eltype.(center)..., eltype.(width)...)
        A = eltype(angle)
        new{S,D,V,C,A,Da}(C.(center), C.(width), angle, value)
    end
end


"""
    Object2d = Object{S,2} where S <: AbstractShape
For 2D objects
"""
const Object2d = Object{S,2} where S

"""
    Object3d = Object{S,3} where S <: AbstractShape
For 3D objects
"""
const Object3d = Object{S,3} where S


# constructors


"""
    Object(shape, center=(0,…), width=(1,…), angle=(0,…), value=1)
    Object(shape ; center, width=(1,…), angle=(0,…), value=1)
General outer `Object` constructor from tuples,
as either positional arguments or named keyword arguments.
"""
function Object(
    shape::AbstractShape{D},
    _center::NTuple{D,RealU} = _tuple(0, D),
    _width::NTuple{D,RealU} = _tuple(1, D),
    _angle::Union{RealU, NTuple{Da,RealU}} where Da = _tuple(0, D-1),
    _value::Number = 1f0 ;
    center::NTuple{D,RealU} = _center,
    width::NTuple{D,RealU} = _width,
    angle::Union{RealU, NTuple{Da,RealU}} where Da = _angle,
    value::Number = _value,
) where {D}
    Object{typeof(shape)}(center, width, angle, value)
end


"""
    Object(shape ; cx, cy, wx=1, wy=wx, ϕ=0, value=1)
    Object(shape ; [6-vector])
2D object constructor from values (without tuples).
"""
function Object(
    shape::AbstractShape{2},
    cx::RealU,
    cy::RealU,
    wx::RealU = oneunit(cx),
    wy::RealU = wx,
    ϕ::RealU = 0,
    value::Number = 1,
) where {C <: RealU}
    Object(shape, (cx, cy), (wx, wy), (ϕ,), value)
end

function Object(shape::AbstractShape{2}, p::AbstractVector ; kwargs...)
    length(p) == 6 || throw(ArgumentError("2D object needs 6 parameters"))
    return Object(shape, p...; kwargs...)
end


"""
    Object(shape ; cx, cy, cz, wx=1, wy=wx, wz=wx, ϕ=0, θ=0, value=1)
    Object(shape ; [9-vector])
3D object constructor from values (without tuples).
"""
function Object(
    shape::AbstractShape{3},
    cx::RealU,
    cy::RealU,
    cz::RealU,
    wx::RealU = oneunit(cx),
    wy::RealU = wx,
    wz::RealU = wx,
    ϕ::RealU = 0,
    θ::RealU = 0,
    value::Number = 1,
) where {C <: RealU}
    Object(shape, (cx, cy, cz), (wx, wy, wz), (ϕ, θ), value)
end

function Object(shape::AbstractShape{3}, p::AbstractVector ; kwargs...)
    length(p) == 9 || throw(ArgumentError("3D object needs 9 parameters"))
    return Object(shape, p...; kwargs...)
end


# methods for objects

"""
    Object(ob::Object ; center, width, angle, value)
Make a copy of Object `ob`, optionally modifying some values.
"""
function Object(ob::Object{S,D} ;
    center::NTuple{D,RealU} = ob.center,
    width::NTuple{D,RealU} = ob.width,
    angle::NTuple{Da,RealU} = ob.angle,
    value::Number = ob.value,
) where {S, D, Da}
    Da == D-1 || throw(ArgumentError("Da=$Da != D-1, where D=$D"))
    Object(S(), center, width, angle, value)
end


Base.eltype(::Object{S,D,V}) where {S,D,V} = V
Base.ndims(::Object{S,D}) where {S,D} = D
Base.ndims(::AbstractShape{D}) where D = D


"""
    show(io::IO, ::MIME"text/plain", ob::Object)
"""
function Base.show(io::IO, ::MIME"text/plain", ob::Object{S,D}) where {S,D}
    println(io, typeof(ob), " (S, D, V, ...)")
    for f in (:center, :width, :angle, :value)
        p = getproperty(ob, f)
        t = typeof(p)
        t = t == NTuple{D,eltype(t)} ? "NTuple{$D,$(eltype(t))}" : "$t"
        println(io, " ", f, "::", t, " ", p)
    end
end


# scale width

"""
    scale(ob::Object, factor::RealU)
    scale(ob::Object, factor::NTuple{D,RealU})
Scale the width(s) by `factor`.
"""
scale(ob::Object{S,D}, factor::NTuple{D,RealU}) where {S,D} =
    Object(S(), ob.center, ob.width .* factor, ob.angle, ob.value)
scale(ob::Object{S,D}, factor::RealU) where {S,D} = scale(ob, _tuple(factor,D))


# scale value

"""
    (*)(ob::Object, x::Number)
    (*)(x::Number, ob::Object)
Scale object `value` by `x`.
"""
Base.:(*)(ob::Object{S}, x::Number) where {S} =
    Object(S(), ob.center, ob.width, ob.angle, ob.value * x)
Base.:(/)(ob::Object, x::Number) = ob * (1 / x)
Base.:(*)(x::Number, ob::Object) = ob * x


# translate

"""
    translate(ob::Object, shift::NTuple{D,RealU})
    translate(ob::Object{S,2}, xshift, yshift)
    translate(ob::Object{S,3}, xshift, yshift, zshift)
Translate the center coordinates of an object by `shift`
"""
translate(ob::Object{S,D}, shift::NTuple{D,RealU}) where {S, D} =
    Object(S(), ob.center .+ shift, ob.width, ob.angle, ob.value)
translate(ob::Object{S,2}, x::RealU, y::RealU) where {S} = translate(ob, (x,y))
translate(ob::Object{S,3}, x::RealU, y::RealU, z::RealU) where {S} = translate(ob, (x,y,z))
