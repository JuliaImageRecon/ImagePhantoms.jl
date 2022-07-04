#=
rect.jl
=#


#using ImagePhantoms #: Object, Object2d

export Rect
export Square
export phantom, radon, spectrum


"""
    Rect <: AbstractShape2
"""
struct Rect <: AbstractShape2 end


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
    end
    return 0*one(T) # zero(T) has units for Unitful :(
end


# methods


"""
    phantom(ob::Object2d{Rect})
Returns function of `(x,y)` for making image.
"""
phantom(ob::Object2d{Rect}) = (x,y) ->
    ob.value * (maximum(abs, coords(ob, x, y)) ≤ 0.5)


"""
    radon_rect(r, ϕ, cx, cy, wx, wy, θ)
Radon transform at `(r,ϕ)` of rectangle (involves a trapezoid).
"""
function radon_rect(r, ϕ, cx, cy, wx, wy, θ)
    (sinϕ, cosϕ) = sincos(ϕ)
    r -= cx * cosϕ + cy * sinϕ # Radon translation property
    (sinϕ, cosϕ) = sincos(ϕ - θ) # Radon rotation property
    xmax = wx * abs(cosϕ)
    ymax = wy * abs(sinϕ)
    lmax = wx * wy / max(xmax, ymax)
    dmax = (xmax + ymax) / 2
    dbreak = abs(xmax - ymax) / 2
    return lmax * trapezoid(r, -dmax, -dbreak, dbreak, dmax)
end

"""
    radon(ob::Object2d{Rect})
Returns function of `(r,ϕ)` for making a sinogram.
"""
radon(ob::Object2d{Rect}) = (r,ϕ) -> ob.value *
    radon_rect(r, ϕ, ob.center..., ob.width..., ob.angle[1])


function spectrum_rect(fx, fy, cx, cy, wx, wy, θ)
    (kx, ky) = rotate2d(fx, fy, θ) # rect is rotated first, then translated
    return cispi(-2 * (fx*cx + fy*cy)) *
           wx * sinc(kx * wx) *
           wy * sinc(ky * wy)
end

"""
    spectrum(ob::Object2d{Rect})
Returns function of ``(f_x,f_y)`` for the spectrum (2D Fourier transform).
"""
spectrum(ob::Object2d{Rect}) = (fx,fy) -> ob.value *
    spectrum_rect(fx, fy, ob.center..., ob.width..., ob.angle[1])
