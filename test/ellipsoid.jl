#=
test/ellipsoid.jl
=#

const DEBUG = false

using ImagePhantoms: Object3d, AbstractShape3, phantom, radon, spectrum
using ImagePhantoms: Ellipsoid, Sphere
import ImagePhantoms as IP
using Unitful: m, unit, °
using FFTW: fftshift, fft
using Test: @test, @testset, @test_throws, @inferred
if DEBUG
    include("helper.jl")
    using MIRTjim: jim, prompt
    using UnitfulRecipes
    using Plots: plot, plot!, scatter, scatter!, gui, default
    default(markerstrokecolor=:auto, markersize=2)
end

(shape, shape3) = (Ellipsoid, Sphere)

macro isob3(ex) # @isob macro to streamline tests
    :(@test $(esc(ex)) isa Object3d{shape})
end


@testset "construct" begin
    @test shape <: AbstractShape3

    # constructors
    @isob3 @inferred Object(shape(), (1,2,3), (4,5,6), (π, π/4), 5.0f0)
    @isob3 @inferred Object(shape(), (1,2,3), (4,5,6), (0, 0), 5.0f0)
    @isob3 @inferred Object(shape(), center=(1,2,3))
    @isob3 @inferred shape((1,2.,3), (4,5//1,6), (π, π/4), 5.0f0)
    @isob3 @inferred shape(1, 2., 3, 4//1, 5, 6., π, π/4, 5.0f0)
    @isob3 @NOTinferred shape(Number[1, 2., 3, 4//1, 5, 6., π, π/4, 5.0f0])

    # spheres
    @isob3 @inferred shape(1, 5.0f0)
    @isob3 @inferred shape3(1, 5.0f0)
    @isob3 @inferred shape3(1, 2, 3, 4., 5.0f0)
    @isob3 @inferred shape3((1, 2, 3), 4., 5.0f0)
    @isob3 @NOTinferred shape3(Number[1, 2, 3, 4., 5.0f0])
end


@testset "operations" begin
    # basic methods

    ob = @inferred shape((1,2.,3), (4,5//1,6), (π, π/4), 5.0f0)

    @isob3 @NOTinferred IP.rotate(ob, π)

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
    x = LinRange(-1,1,32)*5
    y = LinRange(-1,1,31)*5
    z = LinRange(-1,1,30)*5
    ob = @inferred shape((1, 2., 3), (4//1, 5, 6), (π/6, 0), 5.0f0)

    show(devnull, ob)
    @test (@inferred eltype(ob)) == Float32

    fun = @inferred phantom(ob)
    @test fun isa Function
    @test fun(ob.center...) == ob.value
    @test fun((ob.center .+ 2 .* ob.width)...) == 0

    img = @inferred phantom(x, y, z, [ob])
    @test img isa Array{<:Real, 3}
    @test size(img) == length.((x, y, z))

    fun = @inferred radon(ob)
    @test fun isa Function
    fun(0,0,0,0) # todo

    fun = @inferred spectrum(ob)
    @test fun isa Function
    @test fun(0,0,0) ≈ 4/3 * π * prod(ob.width) * ob.value
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
    img = phantom(x, y, z, [ob])

    zscale = 1 / (4/3 * π * prod(width) * ob.value) # normalize spectra
    fx = (-L÷2:L÷2-1) / L / dx
    fy = (-M÷2:M÷2-1) / M / dy
    fz = (-N÷2:N÷2-1) / N / dz
    X = myfft(img) * (dx * dy * dz * zscale)
    kspace = spectrum(fx, fy, fz, [ob]) * zscale
    @test spectrum(ob)(0/m, 0/m, 0/m) * zscale ≈ 1
    @test maximum(abs, kspace) ≈ 1
    @test kspace[L÷2+1,M÷2+1,N÷2+1] ≈ 1

#= todo: move to docs
if DEBUG
    clim = (-6, 0)
    sp = z -> max(log10(abs(z)/oneunit(abs(z))), -6)
    p1 = jim(x, y, img, "phantom")
    p2 = jim(fx, fy, sp.(X), "log10|DFT|"; clim)
    p3 = jim(fx, fy, sp.(kspace), "log10|Spectrum|"; clim)
    p4 = jim(fx, fy, abs.(kspace - X), "Difference")
    jim(p1, p4, p2, p3)
end
=#

    @test abs(maximum(abs, X) - 1) < 1e-2
    @test maximum(abs, kspace - X) / maximum(abs, kspace) < 2e-2

#=
todo
    # test sinogram with projection-slice theorem

    dr = 0.02m
    nr = 2^10
    r = (-nr÷2:nr÷2-1) * dr
    fr = (-nr÷2:nr÷2-1) / nr / dr
    ϕ = deg2rad.(0:180) # * Unitful.rad # todo round unitful Unitful.°
    sino = @inferred radon(r, ϕ, [ob])

    ia = argmin(abs.(ϕ .- deg2rad(55)))
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
=#

end
