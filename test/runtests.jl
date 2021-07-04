# runtests.jl

using Test: @test, @testset, detect_ambiguities
using ImagePhantoms

include("core.jl")

@testset "ImagePhantoms" begin
    @test isempty(detect_ambiguities(ImagePhantoms))
end
