#=
cuboid.jl
=#

include("helper.jl") # todo
using ImagePhantoms #: Object3d, AbstractShape3
using ImagePhantoms #: Cuboid, Cube
import ImagePhantoms as IP
using Unitful: m, unit
#using MIRTjim: jim, prompt
#using UnitfulRecipes
#using Plots; default(markerstrokecolor=:auto)
using FFTW: fftshift, fft
using Test: @test, @testset, @test_throws, @inferred

(shape, shape3) = (Cuboid, Cube)

macro isob3(ex) # macro to streamline tests
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

    @isob3 @inferred shape(1, 5.0f0)
    @isob3 @inferred shape3(1, 5.0f0)
    @isob3 @inferred shape3(1, 2, 3, 4., 5.0f0)
    @isob3 @inferred shape3((1, 2, 3), 4., 5.0f0)
    @isob3 @NOTinferred shape3(Number[1, 2, 3, 4., 5.0f0])
end


@testset "operations" begin
    # basic methods

    ob = @inferred shape((1,2.,3), (4,5//1,6), (π, π/4), 5.0f0)

    @isob @NOTinferred IP.rotate(ob, π)

    @test IP.rotate(ob, -ob.angle[1]).angle[1] == 0

    @isob @inferred ob * 2//1
    @isob @inferred 2 * ob
    @isob @inferred ob / 2.0f0
    @isob @inferred IP.scale(ob, (2,3,4))
    @isob @inferred IP.scale(ob, 2)
    @isob @inferred IP.translate(ob, (2, 3, 4))
    @isob @inferred IP.translate(ob, 2, 3, 4)
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

#=
    fun = @inferred radon(ob)
    @test fun isa Function
    fun(0,0,0,0)
=#

    fun = @inferred spectrum(ob)
    @test fun isa Function
end


@testset "spectrum" begin
    dx = 0.020m
    dy = 0.022m
    dz = 0.025m
    (L,M,N) = (2^8,2^8+2,2^8+4)
    x = (-L÷2:L÷2-1) * dx
    y = (-M÷2:M÷2-1) * dy
    z = (-N÷2:N÷2-1) * dz
    width = (2m, 8m, 3m)
    ob = shape((4m, 3m, 2m), width, (π/6, 0), 1.0f0)
    img = phantom(x, y, z, [ob])

    zscale = 1 / prod(width) # normalize spectra
    fx = (-L÷2:L÷2-1) / L / dx
    fy = (-M÷2:M÷2-1) / M / dy
    fz = (-N÷2:N÷2-1) / N / dz
    X = myfft(img) * (dx * dy * dz * zscale)
    kspace = spectrum(fx, fy, fz, [ob]) * zscale

#=
    clim = (-6, 0)
    sp = z -> max(log10(abs(z)/oneunit(abs(z))), -6)
    p1 = jim(x, y, img, "phantom")
    p2 = jim(fx, fy, sp.(X), "log10|DFT|"; clim)
    p3 = jim(fx, fy, sp.(kspace), "log10|Spectrum|"; clim)
    p4 = jim(fx, fy, abs.(kspace - X), "Difference")
    jim(p1, p4, p2, p3)
=#
    @test maximum(abs, kspace - X) / maximum(abs, kspace) < 2e-2


#=
todo
    # test sinogram with projection-slice theorem

    dr = 0.02m
    nr = 2^10
    r = (-nr÷2:nr÷2-1) * dr
    fr = (-nr÷2:nr÷2-1) / nr / dr
    ϕ = deg2rad.(0:180) # * Unitful.rad # todo round unitful Unitful.°
    sino = radon(r, ϕ, [ob])

    ia = argmin(abs.(ϕ .- deg2rad(55)))
    slice = sino[:,ia]
    Slice = myfft(slice) * dr
    angle = round(rad2deg(ϕ[ia]), digits=1)

    kx, ky = (fr * cos(ϕ[ia]), fr * sin(ϕ[ia])) # Fourier-slice theorem
    ideal = spectrum(ob).(kx, ky)
    @test maximum(abs, ideal - Slice) / maximum(abs, ideal) < 2e-4
=#

#=
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
    plot(p1, p2, p3, p4)
=#

end
