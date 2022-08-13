# test/disk-phantom.jl

using ImagePhantoms: disk_phantom_params, ellipse, phantom, Object2d
using Test: @test, @testset, @inferred

@testset "disk" begin
    params = @inferred disk_phantom_params( ; rhead = () -> rand(100:105))
    @test params isa Vector{<:Tuple}

    ob = @inferred ellipse(params)
    @test ob isa Vector{<:Object2d}
end
