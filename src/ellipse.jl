#=
ellipse.jl
=#

export Ellipse, ellipse, circle


"""
    Ellipse <: AbstractShape{2}
"""
struct Ellipse <: AbstractShape{2} end


# constructor


"""
    ellipse(cx, cy, rx=1, ry=rx, ϕ=0, value::Number=1)
    ellipse(center::NTuple{2,RealU}, radii::NTuple{2,RealU}=(1,1), ϕ::RealU=0, v=1)
Construct `Object{Ellipse}` from parameters.
"""
ellipse(args... ; kwargs...) = Object(Ellipse(), args...; kwargs...)


# special case: circles

"""
    circle(x, y, r, v=1) (circle of radius `r` centered at `(x,y)`)
    circle((x,y), r=1, v=1) ditto
    circle(r, v=1) centered at origin
Construct circles as special cases of `Ellipse`.
"""
circle(cx::RealU, cy::RealU, r::RealU, v::Number = 1) =
    ellipse(cx, cy, r, r, 0, v)
circle(center::NTuple{2,RealU}, r::RealU = oneunit(center[1]), v::Number = 1) =
    circle(center..., r, v)
circle(r::RealU, v::Number = 1) = circle((zero(r), zero(r)), r, v)


# methods


area1(::Ellipse) = π # area of unit circle

ℓmax1(::Ellipse) = 2 # maximum chord through a unit circle


"""
    phantom1(ob::Object2d{Ellipse}, (x,y))
Evaluate unit circle at `(x,y)`,
for unitless coordinates.
"""
phantom1(ob::Object2d{Ellipse}, xy::NTuple{2,Real}) = (sum(abs2, xy) ≤ 1)


# x-ray transform (line integral) of unit circle
# `r` should be unitless
function xray1(::Ellipse, r::R, ϕ::RealU) where {R <: Real}
    T = promote_type(R, Float32)
    r2 = r^2
    return r2 < 1 ? 2 * sqrt(one(T) - r2) : zero(T)
end


"""
    spectrum1(::Object2d{Ellipse}, (kx,ky))
Spectrum of unit circle at `(kx,ky)`,
for unitless spatial frequency coordinates.
"""
function spectrum1(::Object2d{Ellipse}, kxy::NTuple{2,Real})
    return 4 * jinc(2 * sqrt(sum(abs2, kxy)))
end
