# disk-phantom.jl

using ImagePhantoms: disk_phantom_params, ellipse, phantom, Object2d
using Test: @test, @testset, @test_throws, @inferred

@testset "disk" begin
    params = @inferred disk_phantom_params( ; rhead = () -> rand(100:105))
    @test params isa Matrix{Float32}

    ob = @NOTinferred ellipse(params)
    @test ob isa Vector{<:Object2d}
end
