# runtests.jl

using Test: @test, @testset, detect_ambiguities
using ImagePhantoms

include("helper.jl")

include("core.jl")
include("ellipse.jl")
include("gauss2.jl")
include("rect.jl")
include("triangle.jl")
include("shepplogan.jl")
include("disk-phantom.jl")

@testset "ImagePhantoms" begin
    @test isempty(detect_ambiguities(ImagePhantoms))
end
