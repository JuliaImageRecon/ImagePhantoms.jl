


# tests


@show ob = Object(Triangle(), (1,2), (3,4), (5,))
@show ob = Object(Triangle() ; center=(1,2))
@show ob = Object(Triangle(), (1,2), (3,4), (5,))
@show ob = Object(Triangle(), (1,2), (3,4), 5)
Object(Triangle(), (1,2), (3,4), 5., 6., nothing)
Object(Triangle(), (1,2), (3,4), (5.,), 6., nothing)

ob = Object(Triangle() ; center=(1,2))
@show ob
println("here")
rotate(ob, π)
throw(0)

ob = Object(Triangle() ; center=(1,2))
translate(ob, (2,3))
translate(ob, 2, 3)
scale(ob, (2,3))
scale(ob, 2)
throw(0)


#=
abstract type Shape2d <: AbstractShape end
abstract type Shape3d <: AbstractShape end
=#

struct Ellipse <: AbstractShape end
struct Ellipsoid <: AbstractShape end

Object(Ellipse(), (0,0), (1,2))
Object(Ellipse() ; center=(3,4) )
Object(Ellipsoid(), (0,0,0), (1,2,3))
Object(Ellipsoid() ; center=(2,3,4) )

Object(Ellipse(), (0,0), (1,2), π/2, 7, nothing)
Object(Ellipse(), (0,0), (1,2), π/2)
e2 = Object(Ellipse() ; value=7)
# e2 = Object(Ellipse(), 1, 2, 3, 4, π/2, 7)
e2 * 3
e2 / 3
3 * e2
eltype(e2)

translate(e2, -1, 1)
translate(e2, (-1,1))

scale(e2, (2,2))
scale(e2, 2)

rotate(e2, -π)

using Test

#=
=#

#=
struct tmp{D,T}
    a::Dims{D}
    b::NTuple{Da,T} where Da
end
t = @inferred tmp{2,Int}((1,2), (3,))
t = @inferred tmp{3,Int}((1,2,3), (3,4))
typeof(t)


struct Rect <: ShapeType end
struct Gauss <: ShapeType end

"""
    Ellipse(center::NTuple{2,Real}, radii::NTuple{2,Real}, ϕ::Real)
    Ellipse(x, y, rx, ry, ϕ)
    Ellipse(r) (circle of radius `r`)
    circle(x,y,r) (circle of radius `r` centered at `(x,y)`)
    circle((x,y), r) ditto
    circle(r) centered at origin
Construct `Ellipse` from parameters
"""
function Ellipse(center::NTuple{2,Real}, radii::NTuple{2,Real}, ϕ::Real)
	T = promote_type(eltype.(center)..., eltype.(radii)..., eltype(ϕ))
	Shape2d{Ellipse,T}(T.(center), T.(radii), T(ϕ))
end

Ellipse(x::Real, y::Real, rx::Real, ry::Real, ϕ::Real) =
	Ellipse((x, y), (rx, ry), ϕ)

circle(r::Real) = Ellipse((0, 0), (r, r), 0)
circle(x::Real, y::Real, r::Real) = Ellipse((x, y), (r, r), 0)
circle(center::NTuple{2,Real}, r::Real) = Ellipse(center, (r, r), 0)


"""
    Rect(center::NTuple{2,Real}, width::NTuple{2,Real}, ϕ::Real)
    Rect(x, y, rx, ry, ϕ)
    square(x,y,w) (square of width `w` centered at `(x,y)`)
    square((x,y), w) ditto
    square(w) centered at origin
Construct `Rect` from parameters.
"""
function Rect(center::NTuple{2,Real}, width::NTuple{2,Real}, ϕ::Real)
	T = promote_type(eltype.(center)..., eltype.(width)..., eltype(ϕ))
	Shape2d{Rect,T}(T.(center), T.(width), T(ϕ))
end

Rect(x::Real, y::Real, wx::Real, wy::Real, ϕ::Real) =
	Rect((x, y), (wx, wy), ϕ)

square(w::Real) = Rect((0, 0), (w, w), 0)
square(x::Real, y::Real, w::Real) = Rect((x, y), (w, w), 0)
square(center::NTuple{2,Real}, w::Real) = Rect(center, (w, w), 0)


function Gauss(center::NTuple{2,Real}, width::NTuple{2,Real}, ϕ::Real)
	T = promote_type(eltype.(center)..., eltype.(width)..., eltype(ϕ))
	Shape2d{Gauss,T}(T.(center), T.(width), T(ϕ))
end

Gauss(x::Real, y::Real, wx::Real, wy::Real, ϕ::Real) =
	Gauss((x, y), (wx, wy), ϕ)

Gauss(w::Real) = Gauss((0, 0), (w, w), 0)
Gauss(x::Real, y::Real, w::Real) = Gauss((x, y), (w, w), 0)
Gauss(center::NTuple{2,Real}, w::Real) = Gauss(center, (w, w), 0)



# methods





# convert
"""
   convert(shape, T::Type)
Convert `eltype` of `shape` to `T <: Real`
"""
Base.convert(shape::Shape2d{S,O}, T::Type) where {S <: ShapeType, O <: Real} =
	Shape2d{S,T}(T.(shape.center), T.(shape.width), T.(shape.ϕ))


"""
    coordfun(shape::Shape2d, x, y)
Put coordinates `(x,y)` in canonical axes associated with `shape`.
"""
function coordfun(shape::Shape2d{Rect}, x::AbstractArray, y::AbstractArray)
	size(x) == size(y) || throw("x,y size mismatch")
	x = (x .- shape.center[1]) ./ shape.width[1]
	y = (y .- shape.center[2]) ./ shape.width[2]
	c = cos(shape.ϕ)
	s = sin(shape.ϕ)
	return (x * c + y * s, -x * s + y * c) # todo check +/-
end
coordfun(shape::Shape2d{Rect}, x::Real, y::Real) = coordfun(shape, [x], [y])


"""
    imagefun(shape::Shape2d)
Return a function of `(x,y)` for displaying the shape as an image.
"""
function imagefun(shape::Shape2d{Rect})
	return (x,y) -> imagesquare.(coordfun(shape, x, y)...)
end

imagesquare(x,y) = (abs(x) < 1/2) & (abs(y) < 1/2)
imagecirc(x,y) = abs2(x) + abs2(y) < (1/2)^2
#imagegauss(x,y) = exp(-(abs2(x) + abs2(y)) / fwhmfactor) # todo


# tests
using Test

T = Float64
e = Shape2d{Ellipse,T}((1.,2.), (3., 4.), 5.)
e = Shape2d(Ellipse(), 1, 2, 3, 4., π)
e = Shape2d(Ellipse(), (1, 2), (3, 4.), π)
r = Shape2d(Rect(), (1, 2), (3, 4.), π)
@test rotate(r, -r.ϕ).ϕ == 0
e = Ellipse(1, 2, 3, 4., π)
c0 = circle(3)
r = Rect(1, 2, 3, 4., π)
s = square(1)
c1 = translate(c0, (2,4))
c2 = scale(c1, (2,4))
f = convert(e, BigFloat)
eltype(f) == BigFloat

rxy = imagefun(r)
=#
