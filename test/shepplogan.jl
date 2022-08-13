#=
shepplogan.jl
Note: shepp_logan() is type unstable because of possible over-sampling.
=#

using ImagePhantoms # many
import ImagePhantoms as IP
using Test: @test, @testset, @test_throws, @inferred


@testset "shepp" begin
    for case in (SheppLogan, SheppLoganToft, SheppLoganEmis, SheppLoganBrainWeb)
        @test (@inferred case()) isa IP.EllipsePhantomVersion
        @test (@inferred IP.shepp_logan_values(case())) isa Vector
        @test (@inferred ellipse_parameters(case())) isa Vector{<:Tuple}
        image = @inferred shepp_logan(2^5, case())
        type = case == SheppLoganBrainWeb ? Int : Float32
        @test image isa Matrix{type}
    end
    @test (@inferred ellipse_parameters(SouthPark())) isa Vector{<:Tuple}

    image3 = @inferred shepp_logan(2^5, SheppLoganBrainWeb() ; u = (1,2,3))
    @test image3 isa Matrix{Int}

    image0 = @inferred shepp_logan(2^5, SheppLoganEmis())
    image1 = @inferred shepp_logan(2^5, SheppLoganEmis(); oversample=1)
    image2 = @inferred shepp_logan(2^5, SheppLoganEmis(); oversample=3)
    @test image0 == image2

    ob = @inferred shepp_logan(SheppLoganEmis())
    x = LinRange(-1,1,2^5+1) * 0.5
    y = LinRange(-1,1,2^5) * 0.5
    image = @inferred phantom(x, y, ob)
    @test image isa Matrix
    image = phantom(ob).(x,y')
    @test image isa Matrix

    r = LinRange(-0.5,0.5,2^5)
    ϕ = (0:30:180) * deg2rad(1)
    sino = radon(ob).(r,ϕ')
    @test sino isa Matrix
    sino = radon(r, ϕ, ob)
    @test sino isa Matrix

    kx = LinRange(-1,1,2^5) * 9
    ky = LinRange(-1,1,2^5+1) * 9
    kspace = spectrum(ob).(kx, ky')
    @test kspace isa Matrix
    kspace = spectrum(kx, ky, ob)
    @test kspace isa Matrix

    image = @inferred shepp_logan(30, 40, SouthPark(), fovs=(1,1))
end
