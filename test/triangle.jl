#=
triangle.jl
=#

const DEBUG = false

using ImagePhantoms: Object, AbstractShape2, phantom, radon, spectrum
using ImagePhantoms: Triangle
import ImagePhantoms as IP
using Unitful: m, unit
using FFTW: fftshift, fft
using Test: @test, @testset, @test_throws, @inferred
if DEBUG
    include("helper.jl")
    using MIRTjim: jim, prompt; # jim(:prompt,true)
    using UnitfulRecipes
    using Plots; default(markerstrokecolor=:auto, markersize=2)
end

shape = Triangle

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
    @isob @inferred shape(1, 5.0f0)
end


@testset "operations" begin
    # basic methods

    ob = @inferred shape((1,2.), (3,4//1), π, 5.0f0)

    @isob @NOTinferred IP.rotate(ob, π)

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


@testset "spectrum" begin
    dx = 0.02m
    dy = 0.025m
    (M,N) = (2^10,2^10+2)
    x = (-M÷2:M÷2-1) * dx
    y = (-N÷2:N÷2-1) * dy
    width = (13m, 14m)
    ob = shape((-3m, -7m), width, π/6, 1.0f0)
    img = phantom(x, y, [ob])

    zscale = 1 / (sqrt(3)/4 * prod(width)) # normalize spectra by area
    fx = (-M÷2:M÷2-1) / M / dx
    fy = (-N÷2:N÷2-1) / N / dy
    X = myfft(img) * dx * dy * zscale
    kspace = spectrum(fx, fy, [ob]) * zscale

if DEBUG
    clim = (-6, 0)
    sp = z -> max(log10(abs(z)/oneunit(abs(z))), -6)
    p1 = jim(x, y, img, "phantom")
    p2 = jim(fx, fy, sp.(X), "log10|DFT|"; clim)
    p3 = jim(fx, fy, sp.(kspace), "log10|Spectrum|"; clim)
    p4 = jim(fx, fy, abs.(kspace - X), "Difference")
    jim(p1, p4, p2, p3) #; prompt()
end

#   @show maximum(abs, X)
#   @show maximum(abs, kspace)
    @test abs(maximum(abs, X) - 1) < 1e-2
    @test abs(maximum(abs, kspace) - 1) < 1e-5
#   @show maximum(abs, kspace - X) / maximum(abs, kspace)
    @test maximum(abs, kspace - X) / maximum(abs, kspace) < 2e-2


    # test sinogram with projection-slice theorem

    dr = 0.02m
    nr = 2^12
    r = (-nr÷2:nr÷2-1) * dr
    fr = (-nr÷2:nr÷2-1) / nr / dr
    ϕ = deg2rad.(0:360) # * Unitful.rad # todo round unitful Unitful.°
    sino = radon(r, ϕ, [ob])

    ia = argmin(abs.(ϕ .- deg2rad(35)))
    slice = sino[:,ia]
    Slice = myfft(slice) * dr
    angle = round(rad2deg(ϕ[ia]), digits=1)

    kx, ky = (fr * cos(ϕ[ia]), fr * sin(ϕ[ia])) # Fourier-slice theorem
    ideal = spectrum(ob).(kx, ky)

if DEBUG
    p2 = jim(r, rad2deg.(ϕ), sino; aspect_ratio=:none, title="sinogram")
    jim(p1, p2)
    p3 = plot(r, slice, title="profile at ϕ = $angle", label="")
    p4 = scatter(fr, abs.(Slice), label="abs fft", color=:blue)
    scatter!(fr, real(Slice), label="real fft", color=:green)
    scatter!(fr, imag(Slice), label="imag fft", color=:red,
        xlims=(-1,1).*(1.2/m), title="1D spectra")

    plot!(fr, abs.(ideal), label="abs", color=:blue)
    plot!(fr, real(ideal), label="real", color=:green)
    plot!(fr, imag(ideal), label="imag", color=:red)
    plot(p1, p2, p3, p4); gui()
end

    @test maximum(abs, ideal - Slice) / maximum(abs, ideal) < 2e-4
end
