# core.jl

using ImagePhantoms
using ImagePhantoms: rotate2d, coords # helpers
using Unitful: m

using Test: @test, @testset, @test_throws, @inferred

@testset "methods" begin
    ob = Circle(1)
    show(isinteractive() ? stdout : devnull, ob)
    show(isinteractive() ? stdout : devnull, MIME("text/plain"), ob)

    @test ndims(ob) == 2
    @test eltype(ob) == Int64
    @inferred copymod(ob)
    @inferred copymod(ob; value = 2//3)
    @test copymod(ob) == ob
    @test copymod(ob ; value = 3).value == 3

    @test phantom([ob]) isa Function
    @test phantom([ob])(0,0) == 1
    @test phantom(zeros(2), zeros(3), [ob]) == ones(2,3)
    @test phantom(zeros(2,3), zeros(2,3), [ob]) == ones(2,3)

    @test radon([ob]) isa Function
    @test radon([ob])(0,0) == 2
    @test radon(zeros(2), zeros(3), [ob]) == 2*ones(2,3)
    @test radon(zeros(2,3), zeros(2,3), [ob]) == 2*ones(2,3)

    @test spectrum([ob]) isa Function
    @test spectrum([ob])(0,0) ≈ π
    @test spectrum(zeros(2), zeros(3), [ob]) ≈ π*ones(2,3)
    @test spectrum(zeros(2,3), zeros(2,3), [ob]) ≈ π*ones(2,3)

end

@testset "helpers" begin
    @test collect(@inferred rotate2d(2, 1, π/2)) ≈ [1,-2]
    @test all((@inferred rotate2d((2, 1), π/2)) .≈ (1,-2))
    @test (@inferred coords(Square(3m), 9m, 6m)) == (3, 2)
end
