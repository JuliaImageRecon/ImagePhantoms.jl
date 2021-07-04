# core.jl

using ImagePhantoms
using Unitful: m

using Test: @test, @testset, @test_throws, @inferred

@testset "construct" begin
    ell = Ellipse()
    @test ell isa Ellipse
    ig = ImageGeom()
# todo
end

@testset "methods" begin
#   show(isinteractive() ? stdout : devnull, ig)
#   show(isinteractive() ? stdout : devnull, MIME("text/plain"), ig)
# todo
end
