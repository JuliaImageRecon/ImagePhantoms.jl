#=
triangle.jl
2D triangle shape
=#

export Triangle, triangle

const sqrt3 = sqrt(3)


"""
    Triangle{T} <: AbstractShape{2}
By default an equilateral triangle pointing upward with the "center"
in the middle of the base, for the default parameter = `0.5`.

For parameter `p ∈ (0,1)` the object is a triangle whose base is along the x-axis
going from (-p,0) to (1-p,0)` and with height=sqrt(3)/2.

The methods currently support only the default case `p=0.5`.
"""
#struct Triangle{T} <: AbstractShape{2}
struct Triangle <: AbstractShape{2}
#   param::T # fraction in interval (0,1)
#   function Triangle{T}(param::T = 0.5) where {T <: Real}
    function Triangle(param::T = 0.5) where {T <: Real}
        0 < param < 1 || throw(ArgumentError("param=$param"))
        param == 0.5 || throw("Need param = 0.5")
#       new{T}(param)
        new()
    end
end
#Triangle(param::T = 0.5) where {T <: Real} = Triangle{T}(param)
#Triangle(param::T = 0.5) where {T <: Real} = Triangle(param)


# constructors


"""
    triangle(cx, cy, wx=1, wy=wx, ϕ=0, value::Number=1, param::Real=0.5)
    triangle(center::NTuple{2,RealU}, width::NTuple{2,RealU}=(1,1), ϕ::RealU=0, v=1, param=0.5)
Construct `Object{Triangle}` from parameters.
In the typical case where `param=0.5` and `width[1] == width[2]`,
this is an equilateral triangle with base `width[1]` centered along the x axis.
"""
function triangle(args... ; param::Real = 0.5, kwargs...)
    return Object(Triangle(param), args...; kwargs...)
end


# helper

function _trifun(x, y, param)
#   param == 1/2 || throw("") # no need to check because constructor is 1/2
    return (0 ≤ y ≤ sqrt3/2) && y ≤ sqrt3 * (0.5 - abs(x))
end

"""
    (lower, upper) = _interval(a, b)
Determine the interval `[lower, upper]` corresponding to the set
`{x : b ≤ a x}` where `a` is unitless but `b` and `x` have same units.
"""
function  _interval(a::Real, b::RealU)
    Infu = Inf * oneunit(b)
    return a > 0 ? (b/a, Infu) :
        a < 0 ? (-Infu, b/a) :
        b ≤ zero(b) ? (-Infu, Infu) : # a == 0
        (T = typeof(oneunit(b)/one(a)); (zero(T),zero(T)))
end

"""
    radon_tri(r, sinϕ, cosϕ)

For a line integral at radial position `r` and angle `ϕ`,
the locus of points along the line is
`{(r cos(ϕ), r sin(ϕ)) + τ (-sin(ϕ), cos(ϕ)) : τ ∈ ℝ}`.
This function treats the equilateral triangle
with base [-1/2,1/2], pointing upwards,
as the intersection of three half planes:
* H0 = `{(x,y) : y ≥ 0}`
* H1 = `y ≤ √3 (1/2 - x)`
* H2 = `y ≤ √3 (1/2 + x)`.
Find the `τ` values where the line locus lies in each half planes,
then take the length of the intersection of those three intervals.

For example, for H1 we have
`r sin(ϕ) + τ cos(ϕ) ≤ √3 (1/2 - (r cos(ϕ) - τ sin(ϕ)))`
or equivalently
`r (sin(ϕ) + √3 cos(ϕ)) - √3/2 ≤ τ (√3 sin(ϕ) - cos(ϕ))`
which is a set the form `b1 ≤ a1 τ`, corresponding to some interval `(l1,u1)`.
Similarly for H0 and H2.

This approach might not be the most efficient, but it is simple.

See Peter Aundal Toft, "The Radon transform - theory and implementation", 1996
https://orbit.dtu.dk/en/publications/the-radon-transform-theory-and-implementation
for a different approach to finding the Radon transform of a triangle.
"""
function radon_tri(r, sinϕ, cosϕ)
    (l1, u1) = _interval( sqrt3*sinϕ - cosϕ, r*(sinϕ + sqrt3*cosϕ) - sqrt3/2)
    (l2, u2) = _interval(-sqrt3*sinϕ - cosϕ, r*(sinϕ - sqrt3*cosϕ) - sqrt3/2)
    (l0, u0) = _interval(cosϕ, -r*sinϕ)
    tmp = min(u0, u1, u2) - max(l0, l1, l2)
    return tmp > 0 ? tmp : zero(tmp)
end


# spectrum of a unit-base equilateral triangle
function spectrum_tri(u, v)
    (u == 0 && v == 0) ? sqrt3 / 4 :
    (u == 0) ? 1im / (2π * v) * (cispi(-v*sqrt3) - 1) + 1/(2*sqrt3*abs2(π*v)) *
        (1 - cispi(-v * sqrt3) * (1im * π * v * sqrt3 + 1)) :
    1im * sqrt3/(4π*u) * cispi(-sqrt3/2*v) * (
        cispi(-u/2) * sinc((v*sqrt3 - u)/2) -
        cispi( u/2) * sinc((v*sqrt3 + u)/2) )
end


# methods


area1(::Triangle) = sqrt(3)/4 # area of unit base equilateral

ℓmax1(::Triangle) = 1

ℓmax(ob::Object2d{Triangle}) = sqrt(
    (ob.width[1]/2)^2 +
    (ob.width[2] * sqrt(3) / 2)^2
)


"""
    phantom1(ob::Object2d{Triangle}, (x,y))
Evaluate unit triangle at `(x,y)`,
for unitless coordinates.
"""
phantom1(ob::Object2d{Triangle}, xy::NTuple{2,Real}) =
    _trifun(xy..., 0.5)
#   _trifun(xy..., ob.param)


# x-ray transform (line integral) of unit triangle
# `r` should be unitless
function xray1(::Triangle, r::R, ϕ::RealU) where {R <: Real}
    # valid only for param == 1/2, but that's ok now because constructor is 1/2
    T = promote_type(R, Float32)
    return abs(r) ≥ sqrt(3)/2 ? zero(T) : T(radon_tri(r, sincos(ϕ)...))
end


"""
    spectrum1(ob::Object2d{Triangle}, (kx,ky))
Spectrum of unit triangle at `(kx,ky)`,
for unitless spatial frequency coordinates.
"""
function spectrum1(ob::Object2d{Triangle}, kxy::NTuple{2,Real})
#   ob.param == 1/2 || throw("") # no need to check
    return spectrum_tri(kxy...)
end
