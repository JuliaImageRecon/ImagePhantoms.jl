#=
dirac2.jl
=#

export Dirac2, dirac2
#export Dirac


"""
    Dirac{D} <: AbstractShape{D}
"""
struct Dirac{D} <: AbstractShape{D} end

const Dirac2 = Dirac{2}


# constructor


"""
    dirac2(...)
Construct `Object2d{Dirac2}` from parameters;
here `width` is a scale factor: `δ((x - c)/w) = w * δ(x - c)`.
"""
dirac2(args... ; kwargs...) = Object(Dirac2(), args...; kwargs...)
# dirac(D::Int, args... ; kwargs...) = Object(Dirac{D}(), args...; kwargs...)


# methods


area1(::Dirac2) = 1

ℓmax1(::Dirac) = 1 # max line integral


#=
"""
    phantom1(ob::Object2d{Dirac2}, (x,y))
Evaluate Dirac at `(x,y)`,
for unitless coordinates.
Caution: this method is meaningless
because the Dirac impulse is not a function!
"""
phantom1(ob::Object2d{Dirac2}, xy::NTuple{2,Real}) = (xy == (0,0)) ? 1 : 0
=#


#=
"""
    xray1(::Dirac2, r::Real, ϕ::RealU)
X-ray transform (line integral) of Dirac.
`r` should be unitless
Caution: this method is meaningless
because the Dirac impulse is not a function!
"""
xray1(::Dirac2, r::Real, ϕ::RealU) = (r == 0) ? 1 : 0
=#


"""
    spectrum1(::Object2d{Dirac2}, (kx,ky))
Spectrum of Dirac2 at `(kx,ky)`,
for unitless spatial frequency coordinates.
"""
spectrum1(::Object2d{Dirac2}, kxy::NTuple{2,Real}) = 1
#spectrum1(::Object{Dirac, D, V}, kxy::NTuple{D,Real}) where {D, V} = one(V)
