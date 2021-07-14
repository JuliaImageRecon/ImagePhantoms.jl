#=
ellipse.jl
=#


#using ImagePhantoms #: Object2d

export Ellipse
export Circle


"""
    Ellipse <: AbstractShape2
"""
struct Ellipse <: AbstractShape2 end

# constructors


"""
    Ellipse(cx, cy, rx, ry, ϕ, value::Number)
    Ellipse(center::NTuple{2,RealU}, radii::NTuple{2,RealU}, ϕ::RealU, v)
    Ellipse([6-vector])
    Ellipse(r, v=1) (circle of radius `r`)
Construct `Ellipse` object from parameters
"""
function Ellipse(
    cx::RealU,
    cy::RealU,
    rx::RealU,
    ry::RealU,
    ϕ::RealU = 0,
    value::Number = 1,
)
    (cx, cy, rx, ry) = promote(cx, cy, rx, ry)
    Object(Ellipse(), (cx,cy), (rx,ry), ϕ, value)
end

function Ellipse(
    center::NTuple{2,RealU},
    radii::NTuple{2,RealU},
    ϕ::RealU = 0,
    value::Number = 1,
)
    Ellipse(center..., radii..., ϕ, value)
end

function Ellipse(v::AbstractVector{<:Number})
    length(v) == 6 || throw(ArgumentError("$v wrong length"))
    Ellipse(v...)
end

Ellipse(r::RealU, v::Number = 1) = Ellipse((0,0), (r,r), 0, v)


# circles as a special case

"""
    Circle(x,y,r,v=1) (circle of radius `r` centered at `(x,y)`)
    Circle((x,y), r, v=1) ditto
    Circle([4-vector]) ditto
    Circle(r, v=1) centered at origin
Construct circle objects as special cases of `Ellipse` objects.
"""
Circle(r::RealU, v::Number = 1) = Ellipse(r, v)
Circle(cx::RealU, cy::RealU, r::RealU, v::Number = 1) =
    Ellipse(cx, cy,  r, r, 0, v)
Circle(center::NTuple{2,RealU}, r::RealU, v::Number = 1) =
    Ellipse(center, (r, r), 0, v)

function Circle(v::AbstractVector{<:Number})
    length(v) == 4 || throw(ArgumentError("$v wrong length"))
    Circle(v...)
end


# methods

phantom(ob::Object2d{Ellipse}) =
    (x,y) -> (sum(abs2, coords(ob, x, y)) ≤ 1) * ob.value

"""
    radon_ellipse(r, ϕ, cx, cy, rx, ry, θ)
Radon transform of ellipse at point `(r,ϕ)`.
"""
function radon_ellipse(r, ϕ, cx, cy, rx, ry, θ)
	(sinϕ, cosϕ) = sincos(ϕ)
	(sinθ, cosθ) = sincos(θ)
    # square of projected radius:
    rp2 = (rx * (cosϕ * cosθ + sinϕ * sinθ))^2 +
          (ry * (sinϕ * cosθ - cosϕ * sinθ))^2
    sp = cx * cosϕ + cy * sinϕ # radial shift
    dis2 = abs2(r - sp) # square of distances from center
    return 1 / rp2 * sqrt(max(rp2 - dis2, 0))
end

radon(ob::Object2d{Ellipse}) = (r,ϕ) -> ob.value *
    radon_ellipse(r, ϕ, ob.center..., ob.width..., ob.angle[1])

spectrum(ob::Object2d{Ellipse}) = (fx,fy) -> throw("todo")
