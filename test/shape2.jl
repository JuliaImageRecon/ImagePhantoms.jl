# test/shape2.jl

using ImagePhantoms: Object2d, circle, square
using ImagePhantoms: Object, AbstractShape, phantom, radon, spectrum
using ImagePhantoms: Gauss2, gauss2
using ImagePhantoms: Ellipse, ellipse
using ImagePhantoms: Rect, rect
using ImagePhantoms: Triangle, triangle
using ImagePhantoms: Dirac2, dirac2
import ImagePhantoms as IP
using Unitful: m, °
using ImageGeoms: ImageGeom, axesf
using FFTW: fft, fftshift, ifftshift
using Test: @test, @testset, @inferred

include("helper.jl") # myfft

@testset "circle-square" begin # special constructors
    args = [(1, 5.0f0), (1, 2, 3., 5.0f0), ((1, 2), 3., 5.0f0)]
    shapes = (circle, square)
    for shape in shapes, arg in args
        ob = @inferred shape(arg...)
        @test ob isa Object2d
    end
end


@testset "triangle helper" begin
    for (a,b) in [(1, 1), (-1, 1),(0, -1), (0, 1)]
        @inferred IP._interval(a, b)
        @inferred IP._interval(a, b*1m)
    end
    @test IP._interval(2, 4m) == (2m, Inf*m)
end


@testset "gauss2-fwhm" begin
    fwhm = 10
    ob = @inferred gauss2((0, 0), (fwhm, Inf), 0, 1)
    tmp = @inferred phantom((-1:1)*fwhm/2, [0], [ob])
    @test tmp ≈ [0.5, 1, 0.5]
end


# parameters for testing each shape
# (Shape, shape, lmax, lmax1, tol1, tolk, tolp, swidth)
list = [
 (Dirac2, dirac2, 4, 1, Inf, Inf, Inf, (2m, 8m)),
 (Ellipse, ellipse, 8, 2, 1e-6, 6e-4, 2e-4, (2m, 8m)),
 (Gauss2, gauss2, IP.fwhm2spread(4), IP.fwhm2spread(1), 7e-6, 7e-6, 4e-4, (5m, 2m)),
 (Rect, rect, 5, √2, 2e-4, 6e-3, 3e-5, (2m, 8m)),
 (Triangle, triangle, sqrt((4/2)^2 + (3 * sqrt(3) / 2)^2), 1, 3e-5, 2e-3, 3e-5, (13m, 14m)),
]

macro isob(ex) # macro to streamline tests
    :(@test $(esc(ex)) isa Object2d)
end


# constructors
function test2_construct(Shape, shape)
    @test Shape <: AbstractShape{2}
    @isob @inferred Object(Shape(), (1,2), (3,4), π, 5.0f0)
    @isob @inferred Object(Shape(), (1,2), (3,4), (π,), 5.0f0)
    @isob @inferred Object(Shape(), center=(1,2))
    @isob @inferred shape((1,2.), (3,4//1), π, 5.0f0)
    @isob @inferred shape(1, 2., 3, 4//1, π, 5.0f0)
    @isob @inferred shape(1, 5.0f0)
end


# basic methods
function test2_op(Shape, shape)
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


function test2_radon(Shape, shape, ob)
    @inferred IP.xray1(Shape(), 0.5f0, π/6)
    @inferred IP._xray(Shape(), (0., 0.), (2,2), (π/3,), 0.5f0, π/6)

    @test (@inferred IP.xray1(Shape(), 0f0, π/6)) > 0
    @test (@inferred IP._xray(Shape(), (0., 0.), (2,2), (π/3,), 0.5f0, π/6)) > 0

    fun = @inferred radon([ob])
    @inferred fun(0,0)

    r = range(-1,1,51)*2
    s1 = @inferred radon(r, [0], [ob])
    s2 = @inferred radon(r, 0, [ob])
    @test s1[:] == s2

    @test radon([0], [0], [ob])[1] isa Real
    if shape == gauss2
        @test radon([9], [0], [ob])[1] < 1e-1 # outside
    else
        @test radon([9], [0], [ob])[1] == 0 # outside
    end
end


for (Shape, shape, lmax, lmax1, tol1, tolk, tolp, swidth) in list
    @show shape

    has_xray = hasmethod(IP.xray1, (Shape, Real, Real))
    has_phantom = hasmethod(IP.phantom1, (typeof(shape()), NTuple{2,Real}))

    @testset "construct-$shape" begin
        test2_construct(Shape, shape)
    end

    @testset "operations" begin
        test2_op(Shape, shape)
    end

  @testset "method" begin
    x = range(-1,1,51)*5
    y = range(-1,1,50)*5

    ob = @inferred shape((2, 1.), (4//1, 3), π/6, 5.0f0)

    show(devnull, ob)
    @test (@inferred eltype(ob)) == Float32

    @test (@inferred IP.ℓmax(ob)) ≈ lmax
    @test (@inferred IP.ℓmax1(Shape())) ≈ lmax1

    if has_phantom
        fun = @inferred phantom(ob)
        @test fun(ob.center...) == ob.value
        if shape == gauss2
            @test fun((ob.center .+ 9 .* ob.width)...) < 1e-20
        else
            @test fun((ob.center .+ 2 .* ob.width)...) == 0
        end

        img = @inferred phantom(x, y, [ob])
        @test img isa Matrix{<:Real}
        @test size(img) == length.((x, y))
    end

    if has_xray
        test2_radon(Shape, shape, ob)
    end

    fun = @inferred spectrum(ob)
    @test fun isa Function
end


@testset "infer" begin
    obs = [shape((4m, 3m), (2m, 5m), π/6, 1.0f0),
           shape((4f0m, 3f0m), (2f0m, 5f0m), π/6, 1.0)]
    for ob in obs
        nr, dr = 2^4, 0.02m
        r = (-nr÷2:nr÷2-1) * dr .+ ob.center[1]
        if has_xray
            ϕ = (0:30:360) * deg2rad(1)
            @inferred IP._radon(ob, r[1], ϕ[1])
            sino1 = @inferred radon([r[1]], [ϕ[1]], [ob])
            sino = @inferred radon(r, ϕ, [ob])
            @test sino1[1] == sino[1]
        end
    end
  end


  @testset "spectrum" begin
    (M,N) = (2^10,2^10+5) # odd
    ig = ImageGeom(; dims=(M,N), deltas=(0.02m, 0.025m), offsets=:dsp)
    ob = shape((2m, -3m), swidth, π/6, 1.0f0)
    zscale = 1 / (ob.value * IP.area(ob)) # normalize spectra by area

    kspace = (@inferred spectrum(axesf(ig)..., [ob])) * zscale
    @test maximum(abs, kspace) ≈ 1
    @test kspace[M÷2+1,N÷2+1] ≈ 1

    if has_phantom
        img = @inferred phantom(axes(ig)..., [ob])

        X = myfft(img) * (prod(ig.deltas) * zscale)

        @test abs(maximum(abs, X) - 1) < tol1
        errk = maximum(abs, kspace - X) / maximum(abs, kspace)
        @test errk < tolk
    end

    # test sinogram with Fourier-slice theorem

    dr = 0.02m
    nr = 2^10
    r = (-nr÷2:nr÷2-1) * dr
    fr = (-nr÷2:nr÷2-1) / nr / dr
    ϕ = (0:30:360) * deg2rad(1)

    if has_xray
        sino = @inferred radon(r, ϕ, [ob])

        ia = argmin(abs.(ϕ .- deg2rad(50)))
        slice = sino[:,ia]
        Slice = myfft(slice) * dr
        angle = round(rad2deg(ϕ[ia]), digits=1)

        kx, ky = (fr * cos(ϕ[ia]), fr * sin(ϕ[ia])) # Fourier-slice theorem
        ideal = spectrum(ob).(kx, ky)

        errp = maximum(abs, ideal - Slice) / maximum(abs, ideal)
        @test errp < tolp
    end
  end

end # for shape
