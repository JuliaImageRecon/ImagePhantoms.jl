#=
rect.jl
=#

using ImagePhantoms: Object, Object2d, AbstractShape, phantom, radon, spectrum
using ImagePhantoms: Rect, rect
import ImagePhantoms as IP
using Unitful: m, °
using FFTW: fftshift, fft
using Test: @test, @testset, @test_throws, @inferred

(Shape, shape, lmax, lmax1, X1, errk, errps, swidth) =
 (Rect, rect, 5, √2, 2e-4, 6e-3, 3e-5, (2m, 8m))

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

    @test (@inferred IP.ℓmax(ob)) ≈ lmax
    @test (@inferred IP.ℓmax1(Shape())) ≈ lmax1

    fun = @inferred phantom(ob)
    @test fun isa Function
    @test fun(ob.center...) == ob.value
    @test fun((ob.center .+ 2 .* ob.width)...) == 0

    img = @inferred phantom(x, y, [ob])
    @test img isa Matrix{<:Real}
    @test size(img) == length.((x, y))

    @inferred IP.xray1(Shape(), 0.5f0, π/6)
    @inferred IP._xray(Shape(), (0., 0.), (2,2), (π/3,), 0.5f0, π/6)

    fun = @inferred radon([ob])
    @test fun isa Function
    fun(0,0)

    r = LinRange(-1,1,51)*2
    s1 = @inferred radon(r, [0], [ob])
    s2 = @inferred radon(r, 0, [ob])
    @test s1[:] == s2

    fun = @inferred spectrum(ob)
    @test fun isa Function
end


@testset "infer" begin
    obs = [shape((4m, 3m), (2m, 5m), π/6, 1.0f0),
           shape((4f0m, 3f0m), (2f0m, 5f0m), π/6, 1.0)]
    for ob in obs
        nr, dr = 2^4, 0.02m
        r = (-nr÷2:nr÷2-1) * dr .+ ob.center[1]
        ϕ = (0:30:360) * deg2rad(1)
        @inferred IP._radon(ob, r[1], ϕ[1])
        sino1 = @inferred radon([r[1]], [ϕ[1]], [ob])
        sino = @inferred radon(r, ϕ, [ob])
        @test sino1[1] == sino[1]
    end
end


@testset "spectrum" begin
    dx = 0.02m
    dy = 0.025m
    (M,N) = (2^10,2^10+2)
    x = (-M÷2:M÷2-1) * dx
    y = (-N÷2:N÷2-1) * dy
    ob = shape((2m, -3m), swidth, π/6, 1.0f0)
    img = @inferred phantom(x, y, [ob])

    zscale = 1 / (ob.value * IP.area(ob)) # normalize spectra by area
    fx = (-M÷2:M÷2-1) / M / dx
    fy = (-N÷2:N÷2-1) / N / dy
    X = myfft(img) * dx * dy * zscale
    kspace = @inferred spectrum(fx, fy, [ob]) * zscale

    @test abs(maximum(abs, X) - 1) < X1
    @test abs(maximum(abs, kspace) - 1) < 1e-5
    err = maximum(abs, kspace - X) / maximum(abs, kspace)
    @test err < errk


    # test sinogram with projection-slice theorem

    dr = 0.02m
    nr = 2^10
    r = (-nr÷2:nr÷2-1) * dr
    fr = (-nr÷2:nr÷2-1) / nr / dr
    ϕ = (0:30:360) * deg2rad(1)
    sino = @inferred radon(r, ϕ, [ob])

    ia = argmin(abs.(ϕ .- deg2rad(50)))
    slice = sino[:,ia]
    Slice = myfft(slice) * dr
    angle = round(rad2deg(ϕ[ia]), digits=1)

    kx, ky = (fr * cos(ϕ[ia]), fr * sin(ϕ[ia])) # Fourier-slice theorem
    ideal = spectrum(ob).(kx, ky)

    err = maximum(abs, ideal - Slice) / maximum(abs, ideal)
    @test err < errps
end
