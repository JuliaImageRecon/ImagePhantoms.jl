# runtests.jl

using Test: @test, @testset, detect_ambiguities
using ImagePhantoms

include("helper.jl")

include("core.jl")
include("ellipse.jl")
include("gauss2.jl")
include("rect.jl")
include("shepplogan.jl")

@testset "ImagePhantoms" begin
    @test isempty(detect_ambiguities(ImagePhantoms))
end
