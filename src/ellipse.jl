#=
ellipse.jl
=#


#using ImagePhantoms #: Object, Object2d

export Ellipse
export Circle
export phantom, radon, spectrum


"""
    Ellipse <: AbstractShape2
"""
struct Ellipse <: AbstractShape2 end


# constructors


"""
    Ellipse(cx, cy, rx=1, ry=rx, ϕ=0, value::Number=1)
    Ellipse(center::NTuple{2,RealU}, radii::NTuple{2,RealU}=(1,1), ϕ::RealU=0, v=1)
    Ellipse([6-vector])
Construct `Ellipse` object from parameters.
"""
function Ellipse(
    cx::RealU,
    cy::RealU,
    rx::RealU = oneunit(cx),
    ry::RealU = rx,
    ϕ::RealU = 0,
    value::Number = 1,
)
    (cx, cy, rx, ry) = promote(cx, cy, rx, ry)
    Object(Ellipse(), (cx,cy), (rx,ry), ϕ, value)
end

function Ellipse(
    center::NTuple{2,RealU},
    radii::NTuple{2,RealU} = (1,1) .* oneunit(center[1]),
    ϕ::RealU = 0,
    value::Number = 1,
)
    Ellipse(center..., radii..., ϕ, value)
end

function Ellipse(v::AbstractVector{<:Number})
    length(v) == 6 || throw(ArgumentError("$v wrong length"))
    Ellipse(v...)
end


# circles as a special case

"""
    Circle(x,y,r,v=1) (circle of radius `r` centered at `(x,y)`)
    Circle((x,y), r=1, v=1) ditto
    Circle([4-vector]) ditto
    Circle(r, v=1) centered at origin
Construct circle objects as special cases of `Ellipse` objects.
"""
Circle(cx::RealU, cy::RealU, r::RealU, v::Number = 1) =
    Ellipse(cx, cy,  r, r, 0, v)
Circle(center::NTuple{2,RealU}, r::RealU = oneunit(center[1]), v::Number = 1) =
    Circle(center..., r, v)
Circle(r::RealU, v::Number = 1) = Circle((zero(r), zero(r)), r, v)

function Circle(v::AbstractVector{<:Number})
    length(v) == 4 || throw(ArgumentError("$v wrong length"))
    Circle(v...)
end


# helper


# methods


"""
    phantom(ob::Object2d{Ellipse})
Returns function of `(x,y)` for making image.
"""
phantom(ob::Object2d{Ellipse}) =
    (x,y) -> (sum(abs2, coords(ob, x, y)) ≤ 1) * ob.value


"""
    radon_ellipse(r, ϕ, cx, cy, rx, ry, θ)
Radon transform of ellipse at point `(r,ϕ)`.
"""
function radon_ellipse(r, ϕ, cx, cy, rx, ry, θ)
    (sinϕ, cosϕ) = sincos(ϕ)
    r -= cx * cosϕ + cy * sinϕ # Radon translation property
    (sinϕ, cosϕ) = sincos(ϕ - θ) # Radon rotation property
    rp2 = abs2(rx * cosϕ) + abs2(ry * sinϕ) # square of projected radius
    return 2rx*ry / rp2 * sqrt(max(rp2 - abs2(r), 0*oneunit(r)^2))
end


"""
    radon(ob::Object2d{Ellipse})
Returns function of `(r,ϕ)` for making a sinogram.
"""
radon(ob::Object2d{Ellipse}) = (r,ϕ) -> ob.value *
    radon_ellipse(r, ϕ, ob.center..., ob.width..., ob.angle[1])


function spectrum_ellipse(fx, fy, cx, cy, rx, ry, θ)
    (kx, ky) = rotate2d(fx, fy, θ) # rect is rotated first, then translated
    return cispi(-2 * (fx*cx + fy*cy)) * # width=diameter=2*radius
           2rx * 2ry * jinc(2sqrt(abs2(kx * rx) + abs2(ky * ry)))
end

"""
    spectrum(ob::Object2d{Ellipse})
Returns function of ``(f_x,f_y)`` for the spectrum (2D Fourier transform).
"""
spectrum(ob::Object2d{Ellipse}) = (fx,fy) -> ob.value *
    spectrum_ellipse(fx, fy, ob.center..., ob.width..., ob.angle[1])
