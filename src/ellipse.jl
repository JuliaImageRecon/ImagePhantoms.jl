#=
ellipse.jl
=#


#using ImagePhantoms #: Object, Object2d

export Ellipse
export Circle
export phantom, radon, spectrum


"""
    Ellipse <: AbstractShape{2}
"""
struct Ellipse <: AbstractShape{2} end


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
    phantom1(ob::Object2d{Ellipse}, (x,y))
Evaluate unit circle at `(x,y)`,
for unitless coordinates.
"""
phantom1(ob::Object2d{Ellipse}, xy::NTuple{2,Real}) = (sum(abs2, xy) ≤ 1)


# x-ray transform (line integral) of unit circle
# `r` should be unitless
function xray1(::Ellipse, r::Real, ϕ::RealU)
    T = promote_type(eltype(r), Float32)
    r2 = r^2
    return r2 < 1 ? 2 * sqrt(one(T) - r2) : zero(T)
end


"""
    spectrum1(ob::Object2d{Ellipse}, (kx,ky))
Spectrum of unit circle at `(kx,ky)`,
for unitless spatial frequency coordinates.
"""
function spectrum1(ob::Object2d{Ellipse}, kxy::NTuple{2,Real})
    return 4 * jinc(2 * sqrt(sum(abs2, kxy)))
end
