# test/shepplogan.jl

using ImagePhantoms # many
import ImagePhantoms as IP
using Unitful: cm
using Test: @test, @testset, @inferred


@testset "shepp" begin
    (M, N) = (2^5, 2^5+1)
    for case in (SheppLogan, SheppLoganToft, SheppLoganEmis, SheppLoganBrainWeb)
        @test (@inferred case()) isa IP.EllipsePhantomVersion
        @test (@inferred IP.shepp_logan_values(case())) isa Vector
        @test (@inferred ellipse_parameters(case())) isa Vector{<:Tuple}
        image = @inferred shepp_logan(M, case())
        type = case == SheppLoganBrainWeb ? Int : Float32
        @test image isa Matrix{type}
    end
    @test (@inferred ellipse_parameters(SouthPark())) isa Vector{<:Tuple}

    image = @inferred shepp_logan(M, N)
    @test image isa Matrix{Float32}
    @test size(image) == (M, N)

    image = @inferred shepp_logan(M, SheppLoganBrainWeb() ; u = (1,2,3))
    @test image isa Matrix{Int}

    image0 = @inferred shepp_logan(M, SheppLoganEmis())
    image1 = @inferred shepp_logan(M, SheppLoganEmis(); oversample=1)
    image3 = @inferred shepp_logan(M, SheppLoganEmis(); oversample=3)
    @test image0 == image3

    ob = @inferred shepp_logan(SheppLoganEmis())
    x = range(-1,1,M) * 0.5
    y = range(-1,1,N) * 0.5
    image = @inferred phantom(x, y, ob)
    @test image isa Matrix
    image = phantom(ob).(x,y')
    @test image isa Matrix

    r = range(-0.5,0.5,2^5)
    ϕ = (0:30:180) * deg2rad(1)
    sino = radon(ob).(r,ϕ')
    @test sino isa Matrix
    sino = radon(r, ϕ, ob)
    @test sino isa Matrix

    kx = range(-1,1,M) * 9
    ky = range(-1,1,N+1) * 9
    kspace = spectrum(ob).(kx, ky')
    @test kspace isa Matrix
    kspace = spectrum(kx, ky, ob)
    @test kspace isa Matrix

    image = @inferred shepp_logan(M, N, SouthPark(), fovs=(1,1))
end


@testset "shepp3" begin
    @test (@inferred ellipsoid_parameters()) isa Vector{<:Tuple}
    fovs = (8,9,3) .* 1cm
    @inferred ellipsoid_parameters( ; fovs)
    u = (1cm,1,1/cm)
    @inferred ellipsoid_parameters( ; u)
    u = (1,1,1/cm)
    param = @inferred ellipsoid_parameters( ; fovs, u)
    ob = @inferred ellipsoid(param)
    @test ob isa Vector{<:Object3d}
end
