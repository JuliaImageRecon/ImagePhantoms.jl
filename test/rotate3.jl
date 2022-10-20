# test/rotate3.jl

using ImagePhantoms: RealU, Rxyz_mul, Rxyz_inv
using LinearAlgebra: I
using Test: @testset, @test, @inferred


Rx(a) = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)]
Ry(a) = [cos(a) 0 sin(a); 0 1 0; -sin(a) 0 cos(a)]
Rz(a) = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1]

@testset "rotate" begin
    a = randn()
    @test Rx(a)' ≈ Rx(-a)
    @test Rx(a)' ≈ inv(Rx(a))
    @test Ry(a)' ≈ Ry(-a)
    @test Rz(a)' ≈ inv(Rz(a))
    @test Rz(a)' ≈ Rz(-a)
    @test Ry(a)' ≈ inv(Ry(a))
end

function Rxyz_mat(ϕ::RealU, θ::RealU, ψ::RealU)
    return hcat([collect(Rxyz_mul(c..., ϕ, θ, ψ)) for c in eachcol(I(3))]...)
end

function Rxyz_inv_mat(ϕ::RealU, θ::RealU, ψ::RealU)
    return hcat([collect(Rxyz_inv(c..., ϕ, θ, ψ)) for c in eachcol(I(3))]...)
end


@testset "rotate3" begin
    a, b, c = tuple(rand(3)...)
    x, y, z = tuple(randn(Float32, 3)...)

    @test (@inferred Rxyz_mul((x, y, z), a, b, c)) isa NTuple{3, Real}
    @test (@inferred Rxyz_inv((x, y, z), a, b, c)) isa NTuple{3, Real}

    R = Rxyz_mat(a, b, c)
    @test R' ≈ inv(R)

    Ri = Rxyz_inv_mat(a, b, c)
    @test Ri' ≈ inv(Ri)

    @test Ri ≈ inv(R)

    R = Rxyz_mat(a, 0, 0)
    @test R ≈ Rz(a)

    R = Rxyz_mat(0, b, 0)
    @test R ≈ Ry(b)

    R = Rxyz_mat(0, 0, c)
    @test R ≈ Rx(c)

    @test all((@inferred IP.Rxyz_inv((2, 1, 3), π/2, 0, 0)) .≈ (1,-2, 3))
end
