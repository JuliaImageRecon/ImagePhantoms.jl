#=
dirac3.jl
=#

export Dirac3, dirac3

const Dirac3 = Dirac{3}


# constructor


"""
    dirac3(...)
Construct `Object3d{Dirac3}` from parameters;
here `width` is a scale factor: `δ((x - c)/w) = w * δ(x - c)`.
"""
dirac3(args... ; kwargs...) = Object(Dirac3(), args...; kwargs...)


# methods


volume1(::Dirac3) = 1

# ℓmax1(::Dirac) = 1 # max line integral


#=
"""
    phantom1(ob::Object3d{Dirac3}, (x,y))
Evaluate Dirac at `(x,y,z)`,
for unitless coordinates.
Caution: this method is meaningless
because the Dirac impulse is not a function!
"""
phantom1(ob::Object3d{Dirac3}, xyz::NTuple{3,Real}) = (xyz == (0,0,0)) ? 1 : 0
=#


#=
"""
    xray1(::Dirac3, r::Real, ϕ::RealU)
X-ray transform (line integral) of Dirac.
`r` should be unitless
Caution: this method is meaningless
because the Dirac impulse is not a function!
"""
xray1(::Dirac3, u::Real, v::Real, ϕ::RealU, θ::RealU) = ((u,v) == (0,0)) ? 1 : 0
=#


"""
    spectrum1(::Object3d{Dirac3}, (kx,ky,kz))
Spectrum of Dirac3 at `(kx,ky,kz)`,
for unitless spatial frequency coordinates.
"""
spectrum1(::Object3d{Dirac3}, kxyz::NTuple{3,Real}) = 1
