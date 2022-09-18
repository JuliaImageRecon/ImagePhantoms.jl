#=
dirac.jl
=#


using ImagePhantoms: Object, AbstractShape

export Dirac, dirac
export spectrum #phantom, radon


"""
    Dirac{D} <: AbstractShape{D}
"""
struct Dirac{D} <: AbstractShape{D} end


# constructor


"""
    dirac2(cx, cy, wx=1, wy=wx, ϕ=0, value::Number=1)
    dirac2(center::NTuple{2,RealU}, width::NTuple{2,RealU}=(1,1), ϕ::RealU=0, v=1)
Construct `Object{Dirac}` from parameters;
here `width` is the full-width.
"""
#dirac2(args... ; kwargs...) = Object(Dirac{2}(), args...; kwargs...)
#dirac3(args... ; kwargs...) = Object(Dirac{3}(), args...; kwargs...)
dirac(D::Int, args... ; kwargs...) = Object(Dirac{D}(), args...; kwargs...)



# helper


# methods


area1(::Dirac) = 1

ℓmax1(::Dirac) = 1 # max line integral


"""
    phantom1(ob::Object{Dirac}, (x,y))
Evaluate Dirac at `(x,y)`,
for unitless coordinates.
"""
# phantom1(ob::Object{Dirac, D, V}, xy::NTuple{D,Real}) = (xy == (0,0)) ? Inf : zero(V)


# x-ray transform (line integral) of Dirac
# `r` should be unitless
# xray1(::Dirac, r::Real, ϕ::RealU) = (r == 0) ? Inf : zero(r)


"""
    spectrum1(::Object{Dirac}, (kx,ky))
Spectrum of Dirac at `(kx,ky)`,
for unitless spatial frequency coordinates.
"""
spectrum1(::Object{Dirac, D, V}, kxy::NTuple{D,Real}) where {D, V} = one(V)
