#=
object2.jl
Utilities for 2D objects
=#

#export rotate, scale, translate # no: to avoid conflicts with Plots
export phantom, radon, spectrum


area(ob::Object2d{S}) where S = area1(S()) * prod(ob.width)

ℓmax(ob::Object2d{S}) where S = ℓmax1(S()) * maximum(ob.width)


# rotate

"""
    rotate(ob::Object2d, θ::RealU)
Rotate a 2D object.
"""
rotate(ob::Object2d{S}, θ::RealU) where S =
    Object(S(), ob.center, ob.width, ob.angle .+ (θ,), ob.value)


"""
    coords(object::Object2d, x::RealU, y::RealU)
Put coordinates `(x,y)` in canonical axes associated with `object`.
"""
function coords(ob::Object2d, x::RealU, y::RealU)
#   xy = rotate2d(x - ob.center[1], y - ob.center[2], ob.angle[1])
    xy = rotate2d(x - ob.center[1], y - ob.center[2], ob.sin[1], ob.cos[1])
    return xy ./ ob.width # unitless
end


# phantom images


"""
    phantom(ob::Object2d)::Function
Return function of `(x,y)`
that user can sample at any locations
to make a 2D phantom image
for a single 2D object.
"""
function phantom(ob::Object2d)
    # apply coordinate transformation and value:
    return (x,y) -> ob.value * phantom1(ob, coords(ob, x, y))
end


"""
    phantom(oa::Array{<:Object2d})::Function
Return function of `(x,y)`
that user can sample at any locations
to make a 2D phantom image
for an `Array` of 2D objects.
"""
function phantom(oa::Array{<:Object2d})
    return (x,y) -> sum(ob -> phantom(ob)(x,y), oa)
end


"""
    phantom(x::Vector, y::Vector, oa::Array{<:Object2d})
Return a 2D digital image of the phantom sampled at grid of `(x,y)` locations.
Returned 2D image size is `length(x) × length(y)`.
"""
function phantom(
    x::AbstractVector,
    y::AbstractVector,
    oa::Array{<:Object2d},
)
    return phantom(oa).(x, y')
end


"""
    phantom(x::Vector, y::Vector, oa::Array{<:Object2d}, oversample::Int; T)
Return a 2D digital image of the phantom sampled
at grid of `(x,y)` locations,
with over-sampling factor `oversample` and element type `T`.
"""
function phantom(
    x::AbstractVector,
    y::AbstractVector,
    oa::Array{<:Object2d},
    oversample::Int ;
    T::Type{<:Number} = typeof(1f0 * oneunit(promote_type(eltype.(oa)...))),
)

    oversample < 1 && throw(ArgumentError("oversample $oversample"))
    dx = x[2] - x[1]
    dy = y[2] - y[1]
    all(≈(dx), diff(x)) || throw("oversample requires uniform x")
    all(≈(dy), diff(y)) || throw("oversample requires uniform y")
    tmp = ((1:oversample) .- (oversample+1)/2) / oversample
    ophantom = ob -> (x,y) ->
        T(sum(phantom(ob).(x .+ dx*tmp, y .+ dy*tmp')) / abs2(oversample))
    return sum(ob -> ophantom(ob).(x,y'), oa)
end


# radon transform

# shift property of the X-ray transform
function xray_shift(
    r::RealU, ϕ::RealU, # projection coordinates
    cx::RealU, cy::RealU, # object center
)
    (sϕ, cϕ) = sincos(ϕ)
    rshift = cx * cϕ + cy * sϕ
    return r - rshift
end

# rotation property
function xray_rotate(
    Δr::RealU, ϕ::RealU, # projection coordinates
    Φazim::RealU,
)
    return (Δr, ϕ - Φazim)
end

# affine scaling property of the X-ray transform
function xray_scale(
    r::RealU, ϕ::RealU, # projection coordinates
    wx::RealU, wy::RealU, # object width
)
    (sϕ, cϕ) = sincos(ϕ)
    e1old = [cϕ, sϕ]
    e2old = [-sϕ, cϕ]
    width = [wx, wy]
    w = 1 / sqrt( sum(abs2, e2old ./ width) ) # would be radius for circle
    ϕp = atan(wy * sϕ, wx * cϕ) # ϕ'
    tmp = (e1old * r) ./ width
    (sϕp, cϕp) = sincos(ϕp)
    e1new = [cϕp, sϕp]
    rp = e1new' * tmp
    return (w, rp, ϕp)
end


# interface to xray1 after applying shift, rotate, scale properties
function _xray(
    type::AbstractShape{2},
    center::Tuple,
    width::Tuple,
    angle::Tuple,
    r::RealU, ϕ::RealU,
)
    Δr = xray_shift(r, ϕ, center...)
    rr, ϕr = xray_rotate(Δr, ϕ, angle[1])
    scale, rp, ϕp, = xray_scale(rr, ϕr, width...)
    return scale * xray1(type, rp, ϕp)
end


# this gateway seems to help type inference
function _radon(ob::Object2d{S,V,C}, r::Tr, ϕ::Tϕ) where {S,V,C,Tr,Tϕ}
    T = typeof(oneunit(C) * oneunit(V) * one(Tr) * one(Tϕ))
    return T(ob.value * _xray(S(), ob.center, ob.width, ob.angle, r, ϕ))
end


"""
    radon(ob::Object2d)::Function
Return function of `(r,ϕ)` that user can sample
at any projection coordinates
to make sinogram views of a 2D object.

The coordinate system used here is such that `ϕ=0` corresponds to
line integrals along the ``y`` axis
for an object ``f(x,y)``.
Then as `ϕ` increases, the line integrals rotate counter-clockwise.
"""
function radon(ob::Object2d)
    return (r::RealU, ϕ::RealU) -> _radon(ob, r, ϕ)
end


"""
    radon(oa::Array{<:Object2d})::Function
Return function of `(r,ϕ)` that user can sample
at any sinogram coordinates
to make a phantom 2D sinogram
for an array of 3D objects.
"""
function radon(oa::Array{<:Object2d})
    return (r,ϕ) -> sum(ob -> radon(ob)(r,ϕ), oa)
end


"""
    radon(r::Vector, ϕ::Vector, oa::Array{<:Object2d})
Return parallel-beam 2D sinogram
sampled at grid of `(r,ϕ)` locations.
Returned array size is `length(r) × length(ϕ)`.
"""
function radon(
    r::AbstractVector,
    ϕ::AbstractVector,
    oa::Array{<:Object2d},
)
    return radon(oa).(r,ϕ')
end


"""
    radon(r::Vector, ϕ::RealU, oa::Array{<:Object2d})
Return parallel-beam projection
sampled at grid of `r` locations for one `ϕ` value.
Returned array size is `length(r)`.
"""
function radon(
    r::AbstractVector,
    ϕ::RealU,
    oa::Array{<:Object2d},
)
    return radon(oa).(r, ϕ)
end


# spectra


# apply rotate, translate, and scale properties of 2D Fourier transform
function _spectrum(ob::Object2d, fx, fy, cx, cy, rx, ry) #, Φ)
    # rotate then translate:
    (kx, ky) = rotate2d(fx, fy, ob.sin[1], ob.cos[1])
    return rx * ry * cispi(-2*(fx*cx + fy*cy)) *
        spectrum1(ob, (rx*kx, ry*ky))
end


"""
    spectrum(ob::Object2d)::Function
Return function of `(fx,fy)` that user can sample
at any spatial frequency locations
to evaluate the spectrum (2D Fourier transform)
of a single 2D object.
Units of spatial frequencies should be the reciprocal
of the units defining the object.
"""
function spectrum(ob::Object2d)
    return (fx,fy) -> ob.value *
        _spectrum(ob, fx, fy, ob.center..., ob.width...) #, ob.angle[1])
end


"""
    spectrum(oa::Array{<:Object2d})::Function
Return function `kspace(fx,fy)` that user can sample
at any spatial frequency locations
to evaluate the spectrum (2D Fourier transform)
for an `Array` of 2D objects.
"""
function spectrum(oa::Array{<:Object2d})
    return (fx,fy) -> sum(ob -> spectrum(ob)(fx,fy), oa)
end


"""
    spectrum(fx:Vector, fy:Vector, oa::Array{<:Object2d})
Return 2D k-space array
sampled at grid of `(fx,fy)` locations.
Returned 2D array size is `length(fx) × length(fy)`.
"""
function spectrum(
    fx::AbstractVector,
    fy::AbstractVector,
    oa::Array{<:Object2d},
)
    return spectrum(oa).(fx,fy')
end


# helper

function rotate2d(x::RealU, y::RealU, sinθ::Real, cosθ::Real)
    return (cosθ * x + sinθ * y, -sinθ * x + cosθ * y)
end

rotate2d(x::RealU, y::RealU, θ::RealU) = rotate2d(x, y, sincos(θ)...)
rotate2d(xy::NTuple{2,RealU}, args...) = rotate2d(xy..., args...)
