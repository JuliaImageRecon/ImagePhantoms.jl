#=
gauss2.jl
=#

using ImagePhantoms: Object, Object2d, AbstractShape, phantom, radon, spectrum
using ImagePhantoms: Gauss2, gauss2
import ImagePhantoms as IP
using Unitful: m, unit, °
using FFTW: fftshift, fft
using Test: @test, @testset, @test_throws, @inferred

(Shape, shape) = (Gauss2, gauss2)

macro isob(ex) # @isob macro to streamline tests
    :(@test $(esc(ex)) isa Object2d{Shape})
end


@testset "construct" begin
    @test Shape <: AbstractShape{2}

    # constructors
    @isob @inferred Object(Shape(), (1,2), (3,4), π, 5.0f0)
    @isob @inferred Object(Shape(), (1,2), (3,4), (π,), 5.0f0)
    @isob @inferred Object(Shape(), center=(1,2))
    @isob @inferred shape((1,2.), (3,4//1), π, 5.0f0)
    @isob @inferred shape(1, 2., 3, 4//1, π, 5.0f0)
    @isob @inferred shape(1, 5.0f0)
end


@testset "operations" begin
    # basic methods

    ob = @inferred shape((1,2.), (3,4//1), π, 5.0f0)

    @isob @inferred IP.rotate(ob, π)

    @test IP.rotate(ob, -ob.angle[1]).angle[1] == 0

    @isob @inferred ob * 2//1
    @isob @inferred 2 * ob
    @isob @inferred ob / 2.0f0
    @isob @inferred IP.scale(ob, (2,3))
    @isob @inferred IP.scale(ob, 2)
    @isob @inferred IP.translate(ob, (2, 3))
    @isob @inferred IP.translate(ob, 2, 3)
end


@testset "method" begin
    x = LinRange(-1,1,51)*5
    y = LinRange(-1,1,50)*5
    ob = @inferred shape((2, 1.), (4//1, 3), π/6, 5.0f0)

    show(devnull, ob)
    @test (@inferred eltype(ob)) == Float32

    @test (@inferred IP.ℓmax(ob)) ≈ IP.fwhm2spread(4)
    @test (@inferred IP.ℓmax1(Shape())) > 0

    fun = @inferred phantom(ob)
    @test fun isa Function
    @test fun(ob.center...) == ob.value
    @test fun((ob.center .+ 9 .* ob.width)...) < 1e-20

    img = @inferred phantom(x, y, [ob])


    fun = @inferred radon(ob)
    @test fun isa Function
    fun(0,0)

    fun = @inferred spectrum(ob)
    @test fun isa Function
end


@testset "fwhm" begin
    fwhm = 10
    ob = @inferred shape((0, 0), (fwhm, Inf), 0, 1)
    tmp = @inferred phantom((-1:1)*fwhm/2, [0], [ob])
    @test tmp ≈ [0.5, 1, 0.5]
end


@testset "spectrum" begin
    dx = 0.02m
    dy = 0.024m
    (M,N) = (1.5*2^10,2^10+2)
    x = (-M÷2:M÷2-1) * dx
    y = (-N÷2:N÷2-1) * dy
    width = (5m, 2m)
    ob = shape((2m, 3m), width, π/6, 1.0f0)
    img = @inferred phantom(x, y, [ob])

    zscale = 1 / (ob.value * IP.area(ob)) # normalize spectra by area
    fx = (-M÷2:M÷2-1) / M / dx
    fy = (-N÷2:N÷2-1) / N / dy
    X = myfft(img) * dx * dy * zscale
    kspace = @inferred spectrum(fx, fy, [ob]) * zscale

    @test abs(maximum(abs, X) - 1) < 1e-6
    @test abs(maximum(abs, kspace) - 1) < 1e-6
    @test maximum(abs, kspace - X) / maximum(abs, kspace) < 1e-6


    # test sinogram with projection-slice theorem

    dr = 0.02m
    nr = 2^10
    r = (-nr÷2:nr÷2-1) * dr
    fr = (-nr÷2:nr÷2-1) / nr / dr
    ϕ = (0:30:360) * deg2rad(1)
    sino = @NOTinferred radon(r, ϕ, [ob])

    ia = argmin(abs.(ϕ .- deg2rad(55)))
    slice = sino[:,ia]
    Slice = myfft(slice) * dr
    angle = round(rad2deg(ϕ[ia]), digits=1)

    kx, ky = (fr * cos(ϕ[ia]), fr * sin(ϕ[ia])) # Fourier-slice theorem
    ideal = spectrum(ob).(kx, ky)

    @test maximum(abs, ideal - Slice) / maximum(abs, ideal) < 4e-4
end
