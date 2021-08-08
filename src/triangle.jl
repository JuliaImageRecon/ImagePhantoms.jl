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
    Triangle(cx, cy, wx, wy, ϕ, value::Number, param::Real)
    Triangle(center::NTuple{2,RealU}, width::NTuple{2,RealU}, ϕ::RealU, v, param)
    Triangle([6-vector] or [7-vector])
Construct `Triangle` object from parameters.
In the typical case where `param=0.5` and `width[1]==width[2]`,
this is an equilateral triangle with base `width[1]` centered along the x axis.
"""
function Triangle(
    cx::RealU,
    cy::RealU,
    wx::RealU,
    wy::RealU,
    ϕ::RealU = 0,
    value::Number = 1,
    param::Real = 0.5,
)
    (cx, cy, wx, wy) = promote(cx, cy, wx, wy)
    Object(Triangle(), (cx,cy), (wx,wy), ϕ, value, param)
end

function Triangle(
    center::NTuple{2,RealU},
    width::NTuple{2,RealU},
    ϕ::RealU = 0,
    value::Number = 1,
    param::Real = 0.5,
)
    Triangle(center..., radii..., ϕ, value, param)
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


# helper

function _trifun(x, y, param)
    param == 1/2 || throw("todo")
    return (0 ≤ y ≤ sqrt3/2) && abs(x) ≤ 0.5 - y / sqrt3
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
function radon_tri(r, ϕ, cx, cy, wx, wy, θ, p)
    p == 1/2 || throw("todo")
    (sinϕ, cosϕ) = sincos(ϕ)
    r -= cx * cosϕ + cy * sinϕ # Radon translation property
    (sinϕ, cosϕ) = sincos(ϕ - θ) # Radon rotation property
#   xmax = wx * abs(cosϕ)
#   ymax = wy * abs(sinϕ)
#   lmax = wx * wy / max(xmax, ymax)
#   dmax = (xmax + ymax) / 2
#   dbreak = abs(xmax - ymax) / 2
    return NaN # todo
end

"""
    radon(ob::Object2d{Triangle})
Returns function of `(r,ϕ)` for making a sinogram.
"""
radon(ob::Object2d{Triangle}) = (r,ϕ) -> ob.value *
    radon_tri(r, ϕ, ob.center..., ob.width..., ob.angle[1], ob.param)

expimp(x) = exp(1im*π * x)

function spectrum_tri(u, v)
    (u == 0 && v == 0) ? sqrt3 / 4 :
    (u == 0) ? (1 - expimp(-v*sqrt3) / (2π * v)) + 2/sqrt3/(2π*v)^2 *
        (1 - expimp(-v * sqrt3) * ((1im * π * v * sqrt3) - 1)) :
    1/((2π)^2 * u) * (
        expimp(-u) * (1 - expimp(-(v*sqrt3 - u))/(v - u/sqrt3)) -
        expimp( u) * (1 - expimp(-(v*sqrt3 + u))/(v + u/sqrt3)) )
end

function spectrum_tri(fx, fy, cx, cy, wx, wy, θ, param)
    param == 1/2 || throw("todo")
    (kx, ky) = rotate2d(fx, fy, θ) # rotate first, then translated
    return wx * exp(-2im*π*fx*cx) *
           wy * exp(-2im*π*fy*cy) * spectrum_tri(kx/wx, ky/wy)
end

"""
    spectrum(ob::Object2d{Triangle})
Returns function of ``(f_x,f_y)`` for the spectrum (2D Fourier transform).
"""
spectrum(ob::Object2d{Triangle}) = (fx,fy) -> ob.value *
    spectrum_tri(fx, fy, ob.center..., ob.width..., ob.angle[1], ob.param)
