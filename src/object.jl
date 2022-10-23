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
    Object{S, D, V, C, A, Da}(center, width, angle, value) <: AbstractObject
        where {S <: AbstractShape, V <: Number, C,A <: RealU}
General container for 2D and 3D objects for defining image phantoms.

* `center::NTuple{D,C}` coordinates of "center" of this object
* `width::NTuple{D,C}` "width" along each axis
  (e.g., FWHM for Gauss, radii for Ellipse)
* `angle::NTuple{D-1,A}` angle of x' axis relative to x axis,
  in radians (or with units)
* `value::V` "intensity" value for this object

# Example

```jldoctest
julia> Object(Ellipse(), (0,0), (1,2), 0.0, 1//2)
Object2d{Ellipse, Rational{Int64}, Int64, Float64, 1, Float64} (S, D, V, ...)
 center::NTuple{2,Int64} (0, 0)
 width::NTuple{2,Int64} (1, 2)
 angle::Tuple{Float64} (0.0,)
 value::Rational{Int64} 1//2
```
"""
struct Object{S, D, V, C, A, Da, T} <: AbstractObject
    "x,y center coordinates"
    center::NTuple{D,C}
    "'width' along x',y' axes (FWHM for Gauss, radii for Ellipse)"
    width::NTuple{D,C}
    "angle of x' axis relative to x axis, in radians (or with units)"
    angle::NTuple{Da,A}
    "'intensity' value for this shape"
    value::V

    sin::NTuple{Da,T} # sin.(angle)
    cos::NTuple{Da,T} # cos.(angle)

    """
        Object{S}(center, width, angle, value)
    Inner constructor for `S <: AbstractShape`.
    The `center` and `width` tuples should have the same units
    (or should both be unitless).
    """
    function Object{S}(
        center::NTuple{D,RealU},
        width::NTuple{D,RealU},
        angle::NTuple{Da,RealU},
        value::V,
    ) where {S <: AbstractShape, D, Da, V <: Number}
        D == ndims(S()) ||
            throw(ArgumentError("D=$D vs ndims(S)=$(ndims(S)) for S=$S"))
        D == 2 == Da + 1 || D == Da == 3 ||
            throw(ArgumentError("Da=$Da does not fit to D=$D"))
        all(width .> zero(eltype(width))) ||
            throw(ArgumentError("widths must be positive"))

        C = promote_type(eltype.(center)..., eltype.(width)...)
        angle = promote((1f0 .* angle)...) # ensure at least Float32
        A = eltype(angle)
        T = eltype(one(A))
        new{S,D,V,C,A,Da,T}(C.(center), C.(width), angle, value,
            T.(sin.(angle)), T.(cos.(angle)))
    end
end

# handle scalar `angle`
function Object{S}(
    center::NTuple{D,RealU},
    width::NTuple{D,RealU},
    angle::RealU,
    value::V,
) where {S <: AbstractShape, D, V <: Number}
    return Object{S}(center, width, (angle,), value)
end


"""
    Object2d{S,V,C} = Object{S,2,V,C} where {S <: AbstractShape, V,C <: Number}
For 2D objects
"""
const Object2d{S,V,C} = Object{S,2,V,C}

"""
    Object3d{S,V,C} = Object{S,3,V,C} where {S <: AbstractShape, V,C <: Number}
For 3D objects
"""
const Object3d{S,V,C} = Object{S,3,V,C}


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
    _angle::Union{RealU, NTuple{Da,RealU}} = _tuple(0, D == 2 ? 1 : 3),
    _value::Number = 1f0 ;
    center::NTuple{D,RealU} = _center,
    width::NTuple{D,RealU} = _width,
    angle::Union{RealU, NTuple{Da,RealU}} = _angle,
    value::Number = _value,
) where {D,Da}
    Object{typeof(shape)}(center, width, angle, value)
end


"""
    Object(shape ; cx, cy, wx=1, wy=wx, ϕ=0, value=1)
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
)
    Object(shape, (cx, cy), (wx, wy), (ϕ,), value)
end


"""
    Object(shape ; cx, cy, cz, wx=1, wy=wx, wz=wx, ϕ=0, θ=0, ψ=0, value=1)
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
    θ::RealU = zero(ϕ),
    ψ::RealU = zero(ϕ),
    value::Number = 1,
)
    Object(shape, (cx, cy, cz), (wx, wy, wz), (ϕ, θ, ψ), value)
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


"""
    radon_type(::Object)
Determine the element type of the Radon transform of an object
(including units if applicable).
Ensures that its precision is at least `Float32`.
"""
function radon_type(::Object{S, D, V, C, A}) where {S, D, V <: Number, C <: RealU, A <: RealU}
    return eltype(oneunit(C) * oneunit(V) * one(A) * 1f0) # at least Float32
end


"""
    phantom(itr, oa::Array{<:Object})
Return phantom values
sampled at locations
returned by generator (or iterator) `itr`.
Returned array size matches `size(itr)`.
"""
function phantom(
    itr,
    oa::Array{<:Object},
)
    fun = phantom(oa)
    return [fun(i...) for i in itr]
end


"""
    radon(itr, oa::Array{<:Object})
Return parallel-beam projections
sampled at locations
returned by generator (or iterator) `itr`.
Returned array size matches `size(itr)`.
"""
function radon(
    itr,
    oa::Array{<:Object},
)
    fun = radon(oa)
    return [fun(i...) for i in itr]
end


"""
    spectrum(itr, oa::Array{<:Object})
Return spectrum of object(s)
sampled at k-space locations
returned by generator (or iterator) `itr`.
Returned array size matches `size(itr)`.
"""
function spectrum(
    itr,
    oa::Array{<:Object},
)
    fun = spectrum(oa)
    return [fun(i...) for i in itr]
end
