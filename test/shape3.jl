# test/shape3.jl

using ImagePhantoms: Object3d, sphere, cube
using ImagePhantoms: Object, AbstractShape, phantom, radon, spectrum
using ImagePhantoms: Gauss3, gauss3
using ImagePhantoms: Ellipsoid, ellipsoid
using ImagePhantoms: Cone, cone
using ImagePhantoms: Cuboid, cuboid
using ImagePhantoms: Cylinder, cylinder
using ImagePhantoms: Dirac3, dirac3
import ImagePhantoms as IP
using Unitful: m, mm, °
using ImageGeoms: ImageGeom, axesf
using LazyGrids: ndgrid
using FFTW: fft, fftshift, ifftshift
using Test: @test, @testset, @inferred


@testset "sphere-cube" begin # special constructors
    args = [(1, 5.0f0), (1, 2, 3, 4., 5.0f0), ((1, 2, 3), 4., 5.0f0)]
    shapes = (sphere, cube)
    for shape in shapes, arg in args
        ob = @inferred shape(arg...)
        @test ob isa Object3d
    end
end


@testset "xray1 cuboid" begin
    @test (@inferred IP.xray1(Cuboid(), 0, 0, 0, 0)) == 1 # e1,e3
    @test (@inferred IP.xray1(Cuboid(), 0, 0, π/2, 0)) == 1 # e2,e3
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


# parameters for testing each shape
# (Shape, shape, lmax, lmax1, tol1, tolk, tolp)
# for width = (4//1, 5, 6)
list = [
 (Dirac3, dirac3, 6, 1, Inf, Inf, Inf)
 (Ellipsoid, ellipsoid, 12, 2, 1e-2, 4e-4, 1e-3)
 (Gauss3, gauss3, IP.fwhm2spread(6), IP.fwhm2spread(1), 1e-2, 5e-4, 1e-5)
 (Cuboid, cuboid, sqrt(4^2 + 5^2 + 6^2), √3, 1e-2, 2e-2, 2e-3)
 (Cylinder, cylinder, sqrt(10^2 + 6^2), √5, 2e-2, 2e-2, 1e-3)
 (Cone, cone, 10, 2, Inf, Inf, Inf)
]

macro isob3(ex) # macro to streamline tests
    :(@test $(esc(ex)) isa Object3d)
end


# constructors
function test3_construct(Shape, shape)
    @test Shape <: AbstractShape{3}
    @isob3 @inferred Object(Shape(), (1,2,3), (4,5,6), (π, π/4, 0), 5.0f0)
    @isob3 @inferred Object(Shape(), (1,2,3), (4,5,6), (0, 0, 0), 5.0f0)
    @isob3 @inferred Object(Shape(), center=(1,2,3))
    @isob3 @inferred shape((1,2.,3), (4,5//1,6), (π, π/4, 0), 5.0f0)
    @isob3 @inferred shape(1, 2., 3, 4//1, 5, 6., π, π/4, 0, 5.0f0)
end


# basic methods
function test3_op(Shape, shape)
    ob = @inferred shape((1,2.,3), (4,5//1,6), (π, π/4, 0), 5.0f0)

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


function test3_radon(Shape, shape, ob)
    @inferred IP.xray1(Shape(), 0.5f0, 0.3, π/6, π/5)
    @inferred IP._xray(Shape(), (0., 0., 0.), (2,2,2), (π/3,0,0), 0.5f0, 0.3, π/6, π/5)

    @test (@inferred IP.xray1(Shape(), 0, 0, 0, 0)) > 0
    @test (@inferred IP.xray1(Shape(), 0f0, 0., 0, 0)) > 0

    fun = @inferred radon([ob])
    @inferred fun(0,0,0,0)

    @test radon([0], [0], [0], [0], [ob])[1] isa Real
    if shape == gauss3
        @test radon([9], [8], [0], [0], [ob])[1] < 3e-4 # outside
    else
        @test radon([9], [3], [0], [0], [ob])[1] == 0 # outside
    end
end


for (Shape, shape, lmax, lmax1, tol1, tolk, tolp) in list
    @show shape

    has_xray = hasmethod(IP.xray1, (Shape, Real, Real, Real, Real))
    has_phantom = hasmethod(IP.phantom1, (typeof(shape()), NTuple{3,Real}))
    has_spectrum = hasmethod(IP.spectrum1, (typeof(shape()), NTuple{3,Real}))

    @testset "construct-$shape" begin
        test3_construct(Shape, shape)
    end

    @testset "operations" begin
        test3_op(Shape, shape)
    end

  @testset "method" begin
    x = range(-1,1,13)*5
    y = range(-1,1,12)*5
    z = range(-1,1,11)*5

    ob = @inferred shape((1, 2., 3), (4//1, 5, 6), (π/6, 0, 0), 5.0f0)

    show(devnull, ob)
    @test (@inferred eltype(ob)) == Float32

    @test (@inferred IP.ℓmax(ob)) ≈ lmax
    @test (@inferred IP.ℓmax1(Shape())) ≈ lmax1

    if has_phantom
        fun = @inferred phantom(ob)
        @test fun(ob.center...) == ob.value
        if shape == gauss3
            @test fun((ob.center .+ 9 .* ob.width)...) < 1e-20
        else
            @test fun((ob.center .+ 2 .* ob.width)...) == 0
        end

        fun = @inferred phantom([ob])
        @test fun(ob.center...) == ob.value
        if shape == gauss3
            @test fun((ob.center .+ 9 .* ob.width)...) < 1e-20
        else
            @test fun((ob.center .+ 2 .* ob.width)...) == 0
        end

        img = @inferred phantom(x, y, z, [ob])
        @test img isa Array{<:Real, 3}
        @test size(img) == length.((x, y, z))

        over = 2
        img = @inferred phantom(x, y, z, [ob], over)
        @test img isa Array{<:Real, 3}
        @test size(img) == length.((x, y, z))
    end

    volume = IP.volume(ob)

    if has_spectrum
        fun = @inferred spectrum([ob])
        @test (@inferred fun(0,0,0)) ≈ ob.value * volume
    end

    if has_xray
        test3_radon(Shape, shape, ob)

        @testset "radon-units" begin
            width = (30m, 40m, 50m)
            center = (8m, 7m, 6m)
            ϕ = π/6
            ob = shape(center, width, (ϕ, 0, 0), 1.0f0)
            @inferred IP._radon(ob, 0m, (1//2)m, 3f0, 4.)
            x1 = radon([center[1]], [center[3]], [0], [0], [ob])[1]
            x2 = radon([center[1]], [center[3]], 0, 0, [ob])[1]
            @test x2 ≈ x1
        end

    end
  end


  @testset "spectrum" begin
    (L,M,N) = (2^7,2^7+2,2^7+3) # odd
    ig = ImageGeom( dims=(L,M,N), deltas=(1.0m, 1.1m, 1.2m), offsets=:dsp)
    width = (30m, 40m, 50m)
    ob = shape((8m, 7m, 6m), width, (π/6, 0, 0), 5.0f0)
    volume = IP.volume(ob)
    zscale = 1 / (ob.value * volume) # normalize spectra

    if has_spectrum
        kspace = (@inferred spectrum(axesf(ig)..., [ob])) * zscale
        @test spectrum(ob)(0/m, 0/m, 0/m) * zscale ≈ 1
        @test maximum(abs, kspace) ≈ 1
        @test kspace[L÷2+1,M÷2+1,N÷2+1] ≈ 1
    end

    if has_phantom
        oversample = 2
        img = @inferred phantom(axes(ig)..., [ob], oversample)

        X = myfft(img) * (prod(ig.deltas) * zscale)

        @test abs(maximum(abs, X) - 1) < tol1

        if has_spectrum
            errk = maximum(abs, kspace - X) / maximum(abs, kspace)
            @test errk < tolk
        end
    end
  end


  @testset "proj-slice" begin # Fourier-slice theorem
    center = (20mm, 10mm, 5mm)
    width = (25mm, 35mm, 15mm)
    angles = (π/6, 0, 0)
    ob = shape(center, width, angles, 1.0f0)

    pg = ImageGeom((2^8,2^7), (0.6mm,1.0mm), (0.5,0.5)) # projection sampling
    ϕ, θ = π/3, π/7

    if has_xray
        proj = @inferred radon(axes(pg)..., ϕ, θ, [ob])

        e1 = (cos(ϕ), sin(ϕ), 0)
        e3 = (sin(ϕ)*sin(θ), -cos(ϕ)*sin(θ), cos(θ))
        fu, fv = ndgrid(axesf(pg)...)
        ff = vec(fu) * [e1...]' + vec(fv) * [e3...]' # fx,fy,fz for Fourier-slice theorem
        spectrum_slice = spectrum(ob).(ff[:,1], ff[:,2], ff[:,3]) / IP.volume(ob)
        spectrum_slice = reshape(spectrum_slice, pg.dims)
        proj_fft = myfft(proj) * prod(pg.deltas) / IP.volume(ob)
        errp = maximum(abs, spectrum_slice - proj_fft) / maximum(abs, spectrum_slice)
        @test errp < tolp
    end
  end

end # for shape
