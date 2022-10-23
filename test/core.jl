# test/core.jl

using ImagePhantoms: Object, Ellipse, circle, square
using ImagePhantoms: phantom, radon, spectrum
import ImagePhantoms as IP
using Unitful: m, °
using Test: @test, @testset, @inferred


@testset "methods" begin
    ob = circle(1)
    show(isinteractive() ? stdout : devnull, ob)
    show(isinteractive() ? stdout : devnull, MIME("text/plain"), ob)

    @test ndims(ob) == 2
    @test eltype(ob) == Int64
    @inferred Object(ob)
    @inferred Object(ob; value = 2//3)
    @test Object(ob) == ob
    @test Object(ob ; value = 3).value == 3

    @test phantom([ob]) isa Function
    @test phantom([ob])(0,0) == 1
    @test phantom(zeros(2), zeros(3), [ob]) == ones(2,3)

    @test radon([ob]) isa Function
    @test radon([ob])(0,0) == 2
    @test radon(zeros(2), zeros(3), [ob]) == 2*ones(2,3)

    @test spectrum([ob]) isa Function
    @test spectrum([ob])(0.,0.) ≈ π
    @test spectrum(zeros(2), zeros(3), [ob]) ≈ π*ones(2,3)

    @inferred IP.radon_type(ob)
end


@testset "helpers" begin
    @test collect(@inferred IP.rotate2d(2, 1, π/2)) ≈ [1,-2]
    @test all((@inferred IP.rotate2d((2, 1), π/2)) .≈ (1,-2))
    @test (@inferred IP.coords(square(3m), 9m, 6m)) == (3, 2)

    @inferred IP.xray_shift(1.0f0, π/3, 3, 4//5)
    @inferred IP.xray_shift(1.0f0m, π/3, 3m, (4//5)m)
    @inferred IP.xray_rotate(1.0f0, π/3, (1//5))
    @inferred IP.xray_rotate(1.0f0m, π/3, (1//5)°)
    @inferred IP.xray_scale(1.0f0, π/3, 3, 4//5)
    @inferred IP.xray_scale(1.0f0m, π/3, 3m, (4//5)m)
    @inferred IP._xray(Ellipse(), (1.0f0, 2), (1//5, 3.), (π/3,), 2.1, π/7)
    @inferred IP._xray(Ellipse(), (1.0f0m, 2m), ((1//5)m, 3.0m), (π/3,), 2.1m, 5°)
end
