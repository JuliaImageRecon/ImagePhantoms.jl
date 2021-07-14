#=
ellipse.jl
=#

using Revise # todo

using ImagePhantoms #: Object2d, AbstractShape2
using ImagePhantoms #: Ellipse, Circle
using Test: @test, @testset, @test_throws, @inferred

shape = Ellipse

macro isob(ex) # @isob macro to streamline tests
    :(@test $(esc(ex)) isa Object2d{shape})
end

@testset "construct" begin
    @test shape <: AbstractShape2

    # constructors

    ob = @inferred Object(shape(), (1,2), (3,4), π, 5.0f0 )
    @isob ob

    ob = @inferred Object(shape(), (1,2), (3,4), (π,), 5.0f0 )
    @isob ob

    ob = @inferred Object(shape(), center=(1,2))
    @isob ob

    ob = @inferred shape( (1,2.), (3,4//1), π, 5.0f0 )
    @isob ob
    ob = @inferred shape(1, 2., 3, 4//1, π, 5.0f0)
    @isob ob
     ob = shape(Number[1, 2., 3, 4//1, π, 5.0f0]) # @inferred fails
    @isob ob

    # circles
     ob = @inferred shape(1, 5.0f0)
    @isob ob

     ob = @inferred Circle(1, 5.0f0)
    @isob ob
     ob = @inferred Circle(1, 2, 3., 5.0f0)
    @isob ob
     ob = @inferred Circle((1, 2), 3., 5.0f0)
    @isob ob
     ob = Circle(Number[1, 2, 3., 5.0f0]) # @inferred fails
    @isob ob
end

@testset "operations" begin
    # basic methods

    ob = @inferred shape( (1,2.), (3,4//1), π, 5.0f0 )

    oo = rotate(ob, π) # @inferred fails
    @isob oo

    @test rotate(ob, -ob.angle[1]).angle[1] == 0

    oo = @inferred ob * 2//1
    @isob oo

    oo = @inferred 2 * ob
    @isob oo

    oo = @inferred ob / 2.0f0
    @isob oo

    oo = @inferred scale(ob, (2,3))
    @isob oo

    oo = @inferred scale(ob, 2)
    @isob oo

    oo = @inferred translate(ob, (2, 3))
    @isob oo

    oo = @inferred translate(ob, 2, 3)
    @isob oo
end

@testset "method" begin
    ob = @inferred shape( (1,2.), (3,4//1), π/2, 5.0f0 )

    show(devnull, ob)
    @test (@inferred eltype(ob)) == Float32

    fun = @inferred phantom(ob)
    @test fun isa Function
    @test fun(ob.center...) == ob.value
    @test fun((ob.center .+ 2 .* ob.width)...) == 0

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
