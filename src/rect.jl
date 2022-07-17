#=
rect.jl
=#


#using ImagePhantoms #: Object, Object2d

export Rect
export Square
export phantom, radon, spectrum


"""
    Rect <: AbstractShape{2}
"""
struct Rect <: AbstractShape{2} end


# constructors


"""
    Rect(cx, cy, wx=1, wy=wx, ϕ=0, value::Number=1)
    Rect(center::NTuple{2,RealU}, width::NTuple{2,RealU}=(1,1), ϕ::RealU=0, v=1)
    Rect([6-vector])
Construct `Rect` object from parameters;
here `width` is the full-width.
"""
function Rect(
    cx::RealU,
    cy::RealU,
    wx::RealU = oneunit(cx),
    wy::RealU = wx,
    ϕ::RealU = 0,
    value::Number = 1,
)
    (cx, cy, wx, wy) = promote(cx, cy, wx, wy)
    Object(Rect(), (cx,cy), (wx,wy), ϕ, value)
end

function Rect(
    center::NTuple{2,RealU},
    width::NTuple{2,RealU} = (1,1) .* oneunit(center[1]),
    ϕ::RealU = 0,
    value::Number = 1,
)
    Rect(center..., width..., ϕ, value)
end

function Rect(v::AbstractVector{<:Number})
    length(v) == 6 || throw(ArgumentError("$v wrong length"))
    Rect(v...)
end


# squares as a special case

"""
    Square(x,y,w,v=1) (square of width `w` centered at `(x,y)`)
    Square((x,y), w=1, v=1) ditto
    Square([4-vector]) ditto
    Square(w, v=1) centered at origin
Construct `Square` objects as special cases of `Rect` objects.
"""
Square(cx::RealU, cy::RealU, w::RealU, v::Number = 1) =
    Rect(cx, cy, w, w, 0, v)
Square(center::NTuple{2,RealU}, w::RealU = oneunit(center[1]), v::Number = 1) =
    Square(center..., w, v)
Square(w::RealU, v::Number = 1) = Square((zero(w),zero(w)), w, v)

function Square(v::AbstractVector{<:Number})
    length(v) == 4 || throw(ArgumentError("$v wrong length"))
    Square(v...)
end


# helper


"""
    trapezoid(t::Real, t1, t2, t3, t4)
Unit-height trapezoid with breakpoints `t1`, `t2`, `t3`, `t4`.
"""
function trapezoid(t::Real, t1::Real, t2::Real, t3::Real, t4::Real)
    (t, t1, t2, t3, t4) = promote(t, t1, t2, t3, t4)
    T = eltype(t)
    if t1 < t < t2
        return (t - t1) / (t2 - t1)
    elseif t2 ≤ t ≤ t3
        return one(T)
    elseif t3 < t < t4
        return (t4 - t) / (t4 - t3)
    end
    return zero(T)
end


# methods


"""
    phantom1(ob::Object2d{Rect}, (x,y))
Evaluate unit square at `(x,y)`,
for unitless coordinates.
"""
phantom1(ob::Object2d{Rect}, xy::NTuple{2,Real}) = (maximum(abs, xy) ≤ 0.5)


# x-ray transform (line integral) of unit square
# `r` should be unitless
function xray1(::Rect, r::Real, ϕ::RealU)
    (sinϕ, cosϕ) = sincos(ϕ)
    xmax = abs(cosϕ)
    ymax = abs(sinϕ)
    ℓmax = 1 / max(xmax, ymax)
    dmax = (xmax + ymax) / 2
    dbreak = abs(xmax - ymax) / 2
    return ℓmax * trapezoid(r, -dmax, -dbreak, dbreak, dmax)
end


"""
    spectrum1(ob::Object2d{Rect}, (kx,ky))
Spectrum of unit square at `(kx,ky)`,
for unitless spatial frequency coordinates.
"""
function spectrum1(ob::Object2d{Rect}, kxy::NTuple{2,Real})
    return prod(sinc, kxy)
end
