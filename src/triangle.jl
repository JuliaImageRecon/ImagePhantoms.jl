#=
triangle.jl
2D triangle shape
=#


#using ImagePhantoms #: Object, Object2d

export Triangle
export phantom, radon, spectrum

const sqrt3 = sqrt(3)


"""
    Triangle{T} <: AbstractShape2
By default an equilateral triangle pointing upward with the "center"
in the middle of the base, for the default parameter = `0.5`.

For parameter `p` the object is a triangle whose base is along the x-axis
going from (-p,0) to (1-p,0)` and with height=sqrt(3)/2.

Most methods currently support only the case `p=0.5`.
"""
struct Triangle <: AbstractShape2 end

#=
    Triangle{T} <: AbstractShape2
The (width and height can be scaled when constructing an object).
struct Triangle{T} <: AbstractShape2
    param::T # fraction in interval (0,1) 
    function Triangle{T}(param::T = 0.5) where {T <: Real}
        0 < param < 1 || throw(ArgumentError("param=$param"))
        new{T}(param)
    end
end
=#

# constructors


"""
    Triangle(cx, cy, wx=1, wy=wx, ϕ=0, value::Number=1, param::Real=0.5)
    Triangle(center::NTuple{2,RealU}, width::NTuple{2,RealU}=(1,1), ϕ::RealU=0, v=1, param=0.5)
    Triangle([6-vector] or [7-vector])
Construct `Triangle` object from parameters.
In the typical case where `param=0.5` and `width[1] == width[2]`,
this is an equilateral triangle with base `width[1]` centered along the x axis.
"""
function Triangle(
    cx::RealU,
    cy::RealU,
    wx::RealU = oneunit(cx),
    wy::RealU = wx,
    ϕ::RealU = 0,
    value::Number = 1,
    param::Real = 0.5,
)
    (cx, cy, wx, wy) = promote(cx, cy, wx, wy)
    Object(Triangle(), (cx,cy), (wx,wy), ϕ, value, param)
end

function Triangle(
    center::NTuple{2,RealU},
    width::NTuple{2,RealU} = (1,1) .* oneunit(center[1]),
    ϕ::RealU = 0,
    value::Number = 1,
    param::Real = 0.5,
)
    Triangle(center..., width..., ϕ, value, param)
end

function Triangle(v::AbstractVector{<:Number})
    length(v) ∈ (6,7) || throw(ArgumentError("$v wrong length"))
    Triangle(v...)
end

#=
"""
    shape = Triangle(param=0.5)
"""
Triangle(param::T=0.5) where {T <: Real} = Triangle{T}(param)
=#


# helpers

function _trifun(x, y, param)
    param == 1/2 || throw("todo")
    return (0 ≤ y ≤ sqrt3/2) && y ≤ sqrt3 * (0.5 - abs(x))
end

"""
    (lower, upper) = _interval(a, b)
Determine the interval `[lower, upper]` corresponding to the set
`{x : b ≤ a x}`
"""
function  _interval(a::Real, b::RealU)
    return a > 0 ? (b/a, Inf) :
        a < 0 ? (-Inf, b/a) :
        b ≤ 0 ? (-Inf, Inf) : # a == 0
        (T = typeof(one(b)/one(a)); (zero(T),zeros(T)))
end

function radon_tri(r, sinϕ, cosϕ)
    (l1, u1) = _interval( sqrt3*sinϕ - cosϕ, r*(sinϕ + sqrt3*cosϕ) - sqrt3/2)
    (l2, u2) = _interval(-sqrt3*sinϕ - cosϕ, r*(sinϕ - sqrt3*cosϕ) - sqrt3/2)
    (l0, u0) = _interval(cosϕ, -r*sinϕ)
    tmp = min(u0, u1, u2) - max(l0, l1, l2)
    return tmp > 0 ? tmp : zero(tmp)
end


# methods


"""
    phantom(ob::Object2d{Triangle})
Returns function of `(x,y)` for making image.
"""
phantom(ob::Object2d{Triangle}) = (x,y) ->
    ob.value * _trifun(coords(ob, x, y)..., ob.param)


"""
    radon_tri(r, ϕ, cx, cy, wx, wy, θ, p)
Radon transform at `(r,ϕ)` of triangle.
"""
function radon_tri(r, ϕ, cx::C, cy::C, wx::C, wy::C, θ, p) where {C <: RealU}
    T = promote_type(eltype(r), C)
    p == 1/2 || throw("todo")
    (sinϕ, cosϕ) = sincos(ϕ)
    r -= cx * cosϕ + cy * sinϕ # Radon translation property
    ϕ -= θ # Radon rotation property
    (sinϕ, cosϕ) = sincos(ϕ)
    # Radon affine scaling property
    denom = sqrt(abs2(cosϕ/wy) + abs2(sinϕ/wx))
    (sinϕ, cosϕ) = sincos(ϕ)
    r /= wx * wy * denom # unitless
    ϕ = atan(sinϕ/wx, cosϕ/wy)
    return abs(r) ≥ sqrt(3)/2 ? zero(T) : T(radon_tri(r, sincos(ϕ)...) / denom)
end


"""
    radon(ob::Object2d{Triangle})
Returns function of `(r,ϕ)` for making a sinogram.
"""
radon(ob::Object2d{Triangle}) = (r,ϕ) -> ob.value *
    radon_tri(r, ϕ, ob.center..., ob.width..., ob.angle[1], ob.param)


expimp(x::Real) = exp(1im*π * x)

function spectrum_tri(u, v)
    (u == 0 && v == 0) ? sqrt3 / 4 :
    (u == 0) ? 1im / (2π * v) * (expimp(-v*sqrt3) - 1) + 1/(2*sqrt3*abs2(π*v)) *
        (1 - expimp(-v * sqrt3) * (1im * π * v * sqrt3 + 1)) :
    im*sqrt3/(4π*u) * expimp(-sqrt3/2*v) * (
        expimp(-u/2) * sinc((v*sqrt3 - u)/2) -
        expimp( u/2) * sinc((v*sqrt3 + u)/2) )
end

function spectrum_tri(fx, fy, cx, cy, wx, wy, θ, param)
    param == 1/2 || throw("todo")
    (kx, ky) = rotate2d(fx, fy, θ) # rotate first, then translate
    return wx * exp(-2im*π*fx*cx) *
           wy * exp(-2im*π*fy*cy) * spectrum_tri(kx*wx, ky*wy)
end

"""
    spectrum(ob::Object2d{Triangle})
Returns function of ``(f_x,f_y)`` for the spectrum (2D Fourier transform).
"""
spectrum(ob::Object2d{Triangle}) = (fx,fy) -> ob.value *
    spectrum_tri(fx, fy, ob.center..., ob.width..., ob.angle[1], ob.param)
