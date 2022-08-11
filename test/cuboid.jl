#=
test/cuboid.jl
=#

using ImagePhantoms: Object3d, AbstractShape, phantom, radon, spectrum
using ImagePhantoms: Object, Cuboid, cuboid
import ImagePhantoms as IP
using Unitful: m, unit, °
using FFTW: fftshift, fft
using Test: @test, @testset, @test_throws, @inferred

(Shape, shape) = (Cuboid, cuboid)

macro isob3(ex) # @isob macro to streamline tests
    :(@test $(esc(ex)) isa Object3d{Shape})
end


@testset "xray1" begin
    @inferred IP.xray1(Cuboid(), 0.5f0, 0.4, π/6, π/7)
    fun = (u, v, ϕ, θ) -> IP.xray1(Cuboid(), u, v, ϕ, θ)
    @test fun(0, 0, 0, 0) ≈ 1
    @test fun(0.51, 0, 0, 0) == 0
    @test fun(0, 0, π/4, 0) ≈ √2 # diagonal of square
    @test fun(0.72, 0, π/4, 0) == 0
    @test fun(0, 0, π/4, atan(1, √2)) ≈ sqrt(3) # long diagonal of cube
    @inferred IP.cube_bounds(0.2, 0.5f0)
    @inferred IP.cube_bounds(0.2, 0.5)
    @test IP.cube_bounds(0.2, 0.5f0) == IP.cube_bounds(0.2, 0.5)
end


@testset "construct" begin
    @test Shape <: AbstractShape{3}

    # constructors
    @isob3 @inferred Object(Shape(), (1,2,3), (4,5,6), (π, π/4), 5.0f0)
    @isob3 @inferred Object(Shape(), (1,2,3), (4,5,6), (0, 0), 5.0f0)
    @isob3 @inferred Object(Shape(), center=(1,2,3))
    @isob3 @inferred shape((1,2.,3), (4,5//1,6), (π, π/4), 5.0f0)
    @isob3 @inferred shape(1, 2., 3, 4//1, 5, 6., π, π/4, 5.0f0)
end


@testset "operations" begin
    # basic methods

    ob = @inferred shape((1,2.,3), (4,5//1,6), (π, π/4), 5.0f0)

    @isob3 @inferred IP.rotate(ob, π)

    @test IP.rotate(ob, -ob.angle[1]).angle[1] == 0

    @isob3 @inferred ob * 2//1
    @isob3 @inferred 2 * ob
    @isob3 @inferred ob / 2.0f0
    @isob3 @inferred IP.scale(ob, (2,3,4))
    @isob3 @inferred IP.scale(ob, 2)
    @isob3 @inferred IP.translate(ob, (2, 3, 4))
    @isob3 @inferred IP.translate(ob, 2, 3, 4)
end


@testset "method" begin
    x = LinRange(-1,1,13)*5
    y = LinRange(-1,1,12)*5
    z = LinRange(-1,1,11)*5
    ob = @inferred shape((1, 2., 3), (4//1, 5, 6), (π/6, 0), 5.0f0)

    show(devnull, ob)
    @test (@inferred eltype(ob)) == Float32

    @test (@inferred IP.ℓmax(ob)) ≈ sqrt(4^2 + 5^2 + 6^2)
    @test (@inferred IP.ℓmax1(Shape())) ≈ √3

    fun = @inferred phantom(ob)
    @test fun isa Function
    @test fun(ob.center...) == ob.value
    @test fun((ob.center .+ 2 .* ob.width)...) == 0

    fun = @inferred phantom([ob])
    @test fun isa Function
    @test fun(ob.center...) == ob.value
    @test fun((ob.center .+ 2 .* ob.width)...) == 0

    img = @inferred phantom(x, y, z, [ob])
    @test img isa Array{<:Real, 3}
    @test size(img) == length.((x, y, z))

    over = 2
    img = @inferred phantom(x, y, z, [ob], over)
    @test img isa Array{<:Real, 3}
    @test size(img) == length.((x, y, z))


    fun = @inferred radon([ob])
    @test fun isa Function
    fun(0,0,0,0) # todo

    @test radon([0], [0], [0], [0], [ob])[1] isa Real # todo

    volume = IP.volume(ob)

    fun = @inferred spectrum([ob])
    @test fun isa Function
    @test fun(0,0,0) ≈ ob.value * volume
end


@testset "radon-units" begin
    width = (30m, 40m, 50m)
    center = (8m, 7m, 6m)
    ϕ = π/6
    ob = shape(center, width, (ϕ, 0), 1.0f0)
    x1 = radon([center[1]], [center[3]], [0], [0], [ob])[1]
    x2 = radon([center[1]], [center[3]], 0, 0, [ob])[1]
    @test x2 ≈ x1
end


@testset "spectrum" begin
    dx = 1.0m
    dy = 1.1m
    dz = 1.2m
    (L,M,N) = (2^7,2^7+2,2^7+4)
    x = (-L÷2:L÷2-1) * dx
    y = (-M÷2:M÷2-1) * dy
    z = (-N÷2:N÷2-1) * dz
    width = (30m, 40m, 50m)
    ob = shape((8m, 7m, 6m), width, (π/6, 0), 5.0f0)
    oversample = 2
    img = phantom(x, y, z, [ob], oversample)

    volume = IP.volume(ob)
    zscale = 1 / (ob.value * volume) # normalize spectra
    fx = (-L÷2:L÷2-1) / L / dx
    fy = (-M÷2:M÷2-1) / M / dy
    fz = (-N÷2:N÷2-1) / N / dz
    X = myfft(img) * (dx * dy * dz * zscale)
    kspace = spectrum(fx, fy, fz, [ob]) * zscale
    @test spectrum(ob)(0/m, 0/m, 0/m) * zscale ≈ 1
    @test maximum(abs, kspace) ≈ 1
    @test kspace[L÷2+1,M÷2+1,N÷2+1] ≈ 1

    @test abs(maximum(abs, X) - 1) < 1e-2
    err = maximum(abs, kspace - X) / maximum(abs, kspace)
    @test err < 2e-2

    # test sinogram with projection-slice theorem

    du,dv = 0.02m, 0.03m
    nu,nv = 2^9, 2^8
    u = (-nu÷2:nu÷2-1) * du
    v = (-nv÷2:nv÷2-1) * dv
    fu = (-nu÷2:nu÷2-1) / nu / du
    fv = (-nv÷2:nv÷2-1) / nv / dv
#   ϕ = (0:6:180) * deg2rad(1)
    ϕ = deg2rad.(0:6:180) # helps coverage of "abs(e2) < eps(T)"
    θ = [π/7]
    sino = @inferred radon(u, v, ϕ, θ, [ob])

# todo projection slice
end
