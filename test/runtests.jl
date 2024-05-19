# runtests.jl

using Test: @test, @testset, detect_ambiguities
using ImagePhantoms: ImagePhantoms

include("aqua.jl")

include("helper.jl")

include("core.jl")
include("rotate3.jl")

include("shape2.jl")
include("shepplogan.jl")
include("disk-phantom.jl")
include("focus-chart.jl")

include("shape3.jl")

include("iter.jl")

@testset "mri-sense" begin
    include("mri-sense.jl")
end

@testset "ImagePhantoms" begin
    @test isempty(detect_ambiguities(ImagePhantoms))
end
