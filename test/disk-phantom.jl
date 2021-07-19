# disk-phantom.jl

using ImagePhantoms: disk_phantom_params, Ellipse, phantom, Object2d
#using MIRTjim: jim
using Test: @test, @testset, @test_throws, @inferred

@testset "disk" begin
    params = @inferred disk_phantom_params( ; rhead = () -> rand(100:105))
    @test params isa Matrix{Float32}

    ob = Ellipse(params)
    @test ob isa Vector{<:Object2d}
#=
    x = LinRange(-1,1,200)*130
    y = LinRange(-1,1,202)*130
    tmp = phantom(x, y, ob)
    jim(tmp)
=#
end
