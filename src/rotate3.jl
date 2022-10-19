#=
rotate3.jl
3d rotation utilities
=#

const RealU = Real


"""
    Rxyz_mul(x::RealU, y::RealU, z::RealU, ϕ::RealU, θ::RealU, ψ::RealU)
    Rxyz_mul(x::RealU, y::RealU, sinϕ, sinθ, sinψ, cosϕ, cosθ, cosψ)
Multiply `Rx(ψ) * Ry(θ) * Rz(ϕ) * [x, y, z]` for 3D rotation,
where `ψ, θ, ϕ` denote rotation about `x, y, z` axes with right-hand rule.
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


#=
function rotate3d(x::RealU, y::RealU, z::RealU, ϕ::RealU, θ::RealU, ψ::RealU)
    sinϕ, cosϕ = sincos(ψ)
    sinθ, cosθ = sincos(θ)
    sinψ, cosψ = sincos(-ϕ)

    return (
        cosθ * cosψ * x + cosθ * -sinψ * y + sinθ * z,

        (sinϕ * sinθ * cosψ + cosϕ * sinψ) * x +
        (sinϕ * sinθ * -sinψ + cosϕ * cosψ) * y +
        -sinϕ * cosθ * z,

        (cosϕ * -sinθ * cosψ + sinϕ * sinψ) * x +
        (cosϕ * sinθ * sinψ + sinϕ * cosψ) * y +
        cosϕ * cosθ * z,
    )
end

rotate3d(xyz::NTuple{3,RealU}, ϕ::RealU, θ::RealU, ψ::RealU) =
    rotate3d(xyz..., ϕ, θ, ψ)
=#



# tests
# todo - move to test/

Rx(a) = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)]
Ry(a) = [cos(a) 0 sin(a); 0 1 0; -sin(a) 0 cos(a)]
Rz(a) = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1]

using Test
using LinearAlgebra: I

a = randn()
@test Rx(a)' ≈ Rx(-a)
@test Rx(a)' ≈ inv(Rx(a))
@test Ry(a)' ≈ Ry(-a)
@test Rz(a)' ≈ inv(Rz(a))
@test Rz(a)' ≈ Rz(-a)
@test Ry(a)' ≈ inv(Ry(a))

function Rxyz_mat(ϕ::RealU, θ::RealU, ψ::RealU)
    return hcat([collect(Rxyz_mul(c..., ϕ, θ, ψ)) for c in eachcol(I(3))]...)
end

function Rxyz_inv_mat(ϕ::RealU, θ::RealU, ψ::RealU)
    return hcat([collect(Rxyz_inv(c..., ϕ, θ, ψ)) for c in eachcol(I(3))]...)
end

a = randn()
b = randn()
c = randn()

R = Rxyz_mat(a, b, c)
@test R' ≈ inv(R)

Ri = Rxyz_inv_mat(a, b, c)
@test Ri ≈ inv(R)

R = Rxyz_mat(a, 0, 0)
@test R ≈ Rz(a)

R = Rxyz_mat(0, b, 0)
@test R ≈ Ry(b)

R = Rxyz_mat(0, 0, c)
@test R ≈ Rx(c)
