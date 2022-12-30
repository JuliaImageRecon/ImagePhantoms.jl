#=
rect.jl
=#


export Rect, rect, square


"""
    Rect <: AbstractShape{2}
"""
struct Rect <: AbstractShape{2} end


# constructor


"""
    rect(cx, cy, wx=1, wy=wx, ϕ=0, value::Number=1)
    rect(center::NTuple{2,RealU}, width::NTuple{2,RealU}=(1,1), ϕ::RealU=0, v=1)
Construct `Object{Rect}` from parameters;
here `width` is the full-width.
"""
rect(args... ; kwargs...) = Object(Rect(), args...; kwargs...)


# special case: squares

"""
    square(x, y, w,v=1) (square of width `w` centered at `(x,y)`)
    square((x,y), w=1, v=1) ditto
    square(w, v=1) centered at origin
Construct squares as special cases of `Rect`.
"""
square(cx::RealU, cy::RealU, w::RealU, v::Number = 1) =
    rect(cx, cy, w, w, 0, v)
square(center::NTuple{2,RealU}, w::RealU = oneunit(center[1]), v::Number = 1) =
    square(center..., w, v)
square(w::RealU, v::Number = 1) = square((zero(w), zero(w)), w, v)


# helper


"""
    trapezoid(t::Real, t1, t2, t3, t4)
Unit-height trapezoid with breakpoints `t1`, `t2`, `t3`, `t4`.
"""
function trapezoid(t::Real, t1::Real, t2::Real, t3::Real, t4::Real)
    (t, t1, t2, t3, t4) = promote(t, t1, t2, t3, t4)
    T = typeof(t)
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


area1(::Rect) = 1 # area of unit square

ℓmax1(::Rect) = √2 # max line integral through unit square

ℓmax(ob::Object2d{Rect}) = sqrt(sum(abs2, ob.width))


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
    spectrum1(::Object2d{Rect}, (kx,ky))
Spectrum of unit square at `(kx,ky)`,
for unitless spatial frequency coordinates.
"""
function spectrum1(::Object2d{Rect}, kxy::NTuple{2,Real})
    return prod(sinc, kxy)
end
