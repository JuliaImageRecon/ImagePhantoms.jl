#=
test/shape2.jl
=#

using ImagePhantoms: Object2d, circle, square
using Test: @test, @testset, @inferred

@testset "circle-square" begin # special constructors
    args = [(1, 5.0f0), (1, 2, 3., 5.0f0), ((1, 2), 3., 5.0f0)]
    shapes = (circle, square)
    for shape in shapes, arg in args
        ob = @inferred shape(arg...)
        @test ob isa Object2d
    end
end


@testset "triangle helper" begin
    for (a,b) in [(1, 1), (-1, 1),(0, -1), (0, 1)]
        @inferred IP._interval(a, b)
        @inferred IP._interval(a, b*1m)
    end
end
