#=
rotate3.jl
3d rotation utilities
=#


"""
    Rxyz_mul(x::RealU, y::RealU, z::RealU, ϕ::RealU, θ::RealU, ψ::RealU)
    Rxyz_mul(x::RealU, y::RealU, sinϕ, sinθ, sinψ, cosϕ, cosθ, cosψ)
Multiply `Rx(ψ) * Ry(θ) * Rz(ϕ) * [x, y, z]` for 3D rotation,
where `ψ, θ, ϕ` denote rotation about `x, y, z` axes with right-hand rule.

This is the reverse order of z-y′-x″ in
https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Euler_angles_(z-y%E2%80%B2-x%E2%80%B3_intrinsic)_%E2%86%92_rotation_matrix.
"""
function Rxyz_mul(x::RealU, y::RealU, z::RealU,
    sinϕ::Real, sinθ::Real, sinψ::Real,
    cosϕ::Real, cosθ::Real, cosψ::Real,
)
    return (
        cosθ * cosϕ * x +
        cosθ * -sinϕ * y +
        sinθ * z,

        (cosψ * sinϕ + sinψ * sinθ * cosϕ) * x +
        (cosψ * cosϕ - sinψ * sinθ * sinϕ) * y +
        (-sinψ * cosθ) * z,

        (sinψ * sinϕ - cosψ * sinθ * cosϕ) * x +
        (sinψ * cosϕ + cosψ * sinθ * sinϕ) * y +
        cosψ * cosθ * z,
    )
end


"""
    Rxyz_inv(x::RealU, y::RealU, z::RealU, ϕ::RealU, θ::RealU, ψ::RealU)
    Rxyz_inv(x::RealU, y::RealU, sinϕ, sinθ, sinψ, cosϕ, cosθ, cosψ)
Multiply `Rz(-ϕ) * Ry(-θ) * Rx(-ψ) * [x, y, z]` for 3D (inverse) rotation.
"""
function Rxyz_inv(x::RealU, y::RealU, z::RealU,
    sinϕ::Real, sinθ::Real, sinψ::Real,
    cosϕ::Real, cosθ::Real, cosψ::Real,
)
    return (
        cosθ * cosϕ * x +
        (cosψ * sinϕ + sinψ * sinθ * cosϕ) * y +
        (sinψ * sinϕ - cosψ * sinθ * cosϕ) * z,

        cosθ * -sinϕ * x +
        (cosψ * cosϕ - sinψ * sinθ * sinϕ) * y +
        (sinψ * cosϕ + cosψ * sinθ * sinϕ) * z,

        sinθ * x +
        (-sinψ * cosθ) * y +
        cosψ * cosθ * z,
    )
end


function Rxyz_mul(x::RealU, y::RealU, z::RealU, ϕ::RealU, θ::RealU, ψ::RealU)
    sinϕ, cosϕ = sincos(ϕ)
    sinθ, cosθ = sincos(θ)
    sinψ, cosψ = sincos(ψ)
    return Rxyz_mul(x, y, z, sinϕ, sinθ, sinψ, cosϕ, cosθ, cosψ)
end


function Rxyz_inv(x::RealU, y::RealU, z::RealU, ϕ::RealU, θ::RealU, ψ::RealU)
    sinϕ, cosϕ = sincos(ϕ)
    sinθ, cosθ = sincos(θ)
    sinψ, cosψ = sincos(ψ)
    return Rxyz_inv(x, y, z, sinϕ, sinθ, sinψ, cosϕ, cosθ, cosψ)
end


Rxyz_mul(xyz::NTuple{3,RealU}, args...) = Rxyz_mul(xyz..., args...)

Rxyz_inv(xyz::NTuple{3,RealU}, args...) = Rxyz_inv(xyz..., args...)
