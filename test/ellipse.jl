#=
ellipse.jl
=#

using Revise # todo

using ImagePhantoms #: Object2d, AbstractShape2
using ImagePhantoms #: Ellipse, Circle
using Test: @test, @testset, @test_throws, @inferred
include("infer.jl")

(shape, shape2) = (Ellipse, Circle)

macro isob(ex) # @isob macro to streamline tests
    :(@test $(esc(ex)) isa Object2d{shape})
end


@testset "construct" begin
    @test shape <: AbstractShape2

    # constructors
    @isob @inferred Object(shape(), (1,2), (3,4), π, 5.0f0)
    @isob @inferred Object(shape(), (1,2), (3,4), (π,), 5.0f0)
    @isob @inferred Object(shape(), center=(1,2))
    @isob @inferred shape((1,2.), (3,4//1), π, 5.0f0)
    @isob @inferred shape(1, 2., 3, 4//1, π, 5.0f0)
    @isob @NOTinferred shape(Number[1, 2., 3, 4//1, π, 5.0f0])

    # circles
    @isob @inferred shape(1, 5.0f0)
    @isob @inferred shape2(1, 5.0f0)
    @isob @inferred shape2(1, 2, 3., 5.0f0)
    @isob @inferred shape2((1, 2), 3., 5.0f0)
    @isob @NOTinferred shape2(Number[1, 2, 3., 5.0f0])
end


@testset "operations" begin
    # basic methods

    ob = @inferred shape((1,2.), (3,4//1), π, 5.0f0)

    @isob @NOTinferred rotate(ob, π)

    @test rotate(ob, -ob.angle[1]).angle[1] == 0

    @isob @inferred ob * 2//1
    @isob @inferred 2 * ob
    @isob @inferred ob / 2.0f0
    @isob @inferred scale(ob, (2,3))
    @isob @inferred scale(ob, 2)
    @isob @inferred translate(ob, (2, 3))
    @isob @inferred translate(ob, 2, 3)
end


@testset "method" begin
	x = LinRange(-1,1,51)*5
    y = LinRange(-1,1,50)*5
    ob = @inferred shape((1,2.), (3,4//1), π/2, 5.0f0)

    show(devnull, ob)
    @test (@inferred eltype(ob)) == Float32

    fun = @inferred phantom(ob)
    @test fun isa Function
    @test fun(ob.center...) == ob.value
    @test fun((ob.center .+ 2 .* ob.width)...) == 0

    img = @inferred phantom(x, y, [ob])
    fun = @inferred radon(ob)
    @test fun isa Function
    fun(0,0)

    fun = @inferred spectrum(ob)
    @test fun isa Function
end

#=
todo: test with projection-slice theorem
    r = LinRange(-1,1,1001) * 8
    ϕ = LinRange(0,π,1003)
    sino = fun.(r, ϕ')
    jim(r, ϕ, sino)

todo: test FFT vs analytical, at least if is interactive
=#
