# test/focus-chart.jl

using ImagePhantoms: focus_chart, Object2d, Triangle, phantom
using Test: @test, @testset, @inferred


@testset "focus-chart" begin
    ob = @inferred focus_chart( ; nspoke=20)
    @test ob isa Vector{<: Object2d{Triangle}}

    x = range(-1,1,2^9) * 1.1
    image = phantom(ob).(x,x')
    @test image isa Matrix
end
