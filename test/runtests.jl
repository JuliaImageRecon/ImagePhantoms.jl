# runtests.jl

using Test: @test, @testset, detect_ambiguities
using ImagePhantoms

include("helper.jl")

include("core.jl")

# 2d:
include("ellipse.jl")
include("gauss2.jl")
include("rect.jl")
include("triangle.jl")
include("shepplogan.jl")
include("disk-phantom.jl")
include("focus-chart.jl")

# 3d:
include("ellipsoid.jl")
include("gauss3.jl")
include("cuboid.jl")
include("cylinder.jl")

@testset "mri-sense" begin
    include("mri-sense.jl")
end

@testset "ImagePhantoms" begin
    @test isempty(detect_ambiguities(ImagePhantoms))
end
