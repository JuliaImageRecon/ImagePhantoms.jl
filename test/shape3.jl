#=
test/shape3.jl
=#

using ImagePhantoms: Object3d, sphere, cube
using Test: @test, @testset, @inferred

@testset "sphere-cube" begin # special constructors
    args = [(1, 5.0f0), (1, 2, 3, 4., 5.0f0), ((1, 2, 3), 4., 5.0f0)]
    shapes = (sphere, cube)
    for shape in shapes, arg in args
        ob = @inferred shape(arg...)
        @test ob isa Object3d
    end
end
