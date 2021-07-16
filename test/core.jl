# core.jl

using ImagePhantoms
using ImagePhantoms: rotate2d, coords # helpers
using Unitful: m

using Test: @test, @testset, @test_throws, @inferred

@testset "construct" begin
    ell = Ellipse()
    @test ell isa Ellipse
#   ig = ImageGeom()
# todo
end

@testset "methods" begin
#   show(isinteractive() ? stdout : devnull, ig)
#   show(isinteractive() ? stdout : devnull, MIME("text/plain"), ig)
# todo
end

@testset "helpers" begin
    @test collect(@inferred rotate2d(2, 1, π/2)) ≈ [1,-2]
    @test (@inferred coords(Square(3m), 9m, 6m)) == (3, 2)
end
