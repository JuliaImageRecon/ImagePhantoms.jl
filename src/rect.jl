#=
rect.jl
=#


#using ImagePhantoms #: Object2d

export Rect
export Square


"""
    Rect <: AbstractShape2
"""
struct Rect <: AbstractShape2 end


# constructors


"""
    Rect(cx, cy, wx, wy, ϕ, value::Number)
    Rect(center::NTuple{2,RealU}, width::NTuple{2,RealU}, ϕ::RealU, v)
    Rect([6-vector])
    Rect(r, v=1) (square of radius `r`)
Construct `Rect` object from parameters
"""
function Rect(
    cx::RealU,
    cy::RealU,
    wx::RealU,
    wy::RealU,
    ϕ::RealU = 0,
    value::Number = 1,
)
    (cx, cy, wx, wy) = promote(cx, cy, wx, wy)
    Object(Rect(), (cx,cy), (wx,wy), ϕ, value)
end

function Rect(
    center::NTuple{2,RealU},
    width::NTuple{2,RealU},
    ϕ::RealU = 0,
    value::Number = 1,
)
    Rect(center..., width..., ϕ, value)
end

function Rect(v::AbstractVector{<:Number})
    length(v) == 6 || throw(ArgumentError("$v wrong length"))
    Rect(v...)
end

Rect(r::RealU, v::Number = 1) = Rect((0,0), (r,r), 0, v)


# squares as a special case

"""
    Square(x,y,w,v=1) (square of width `w` centered at `(x,y)`)
    Square((x,y), w, v=1) ditto
    Square([4-vector]) ditto
    Square(w, v=1) centered at origin
Construct `Square` objects as special cases of `Rect` objects.
"""
Square(w::RealU, v::Number = 1) = Rect(w, v)
Square(cx::RealU, cy::RealU, w::RealU, v::Number = 1) =
    Rect(cx, cy, w, w, 0, v)
Square(center::NTuple{2,RealU}, w::RealU, v::Number = 1) =
    Rect(center, (w, w), 0, v)

function Square(v::AbstractVector{<:Number})
    length(v) == 4 || throw(ArgumentError("$v wrong length"))
    Square(v...)
end


# helpers


"""
    trapezoid(t::RealU, t1, t2, t3, t4)
Unit-height trapezoid with breakpoints `t1`, `t2`, `t3`, `t4`.
"""
function trapezoid(t::RealU, t1::RealU, t2::RealU, t3::RealU, t4::RealU)
	(t, t1, t2, t3, t4) = promote(t, t1, t2, t3, t4)
	T = eltype(t)
	if t1 < t < t2
		return (t - t1)/(t2 - t1)
	elseif t2 <= t <= t3
		return one(T)
	elseif t3 < t < t4
		return (t4 - t)/(t4 - t3)
	else
		return zero(T)
	end
end


# methods


phantom(ob::Object2d{Rect}) = (x,y) ->
    ob.value * (sum(c -> max(abs.(c)), coords(ob, x, y)) ≤ 1)


"""
    radon_rect(r, ϕ, cx, cy, wx, wy, θ)
Radon transform of rectangle at point `(r,ϕ)`.
"""
function radon_rect(r, ϕ, cx, cy, wx, wy, θ)
    (sinϕ, cosϕ) = sincos(ϕ + θ) # Radon rotation property, todo: check sign
    xmax = wx * abs(cos(ϕ))
    ymax = wy * abs(sin(ϕ))
    lmax = wx * wy / maximum(xmax, ymax)
    dmax = (xmax + ymax) / 2
    dbreak = abs(xmax - ymax) / 2
    return trapezoid(r, -dmax, -dbreak, dbreak, dmax)
end

radon(ob::Object2d{Rect}) = (r,ϕ) -> ob.value *
    radon_rect(r, ϕ, ob.center..., ob.width..., ob.angle[1])


function spectrum_rect(fx, fy, cx, cy, wx, wy, θ)
    (kx, ky) = rotate2d(fx, fy, θ)
    return wx * sinc(kx * wx) * exp(-2im*π*kx*cx) *
           wy * sinc(ky * wy) * exp(-2im*π*ky*cy)
end

spectrum(ob::Object2d{Rect}) = (fx,fy) -> ob.value *
    spectrum_rect(fx, fy, ob.center..., ob.width..., ob.angle[1])
