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
include("focus-chart.jl")

@testset "mri-sense" begin
    include("mri-sense.jl")
end

@testset "ImagePhantoms" begin
    @test isempty(detect_ambiguities(ImagePhantoms))
end
