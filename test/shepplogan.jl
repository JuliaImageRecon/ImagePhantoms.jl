#=
shepplogan.jl
=#

using ImagePhantoms # many
using Test: @test, @testset, @test_throws, @inferred


@testset "shepp" begin
    for case in (SheppLogan, SheppLoganToft, SheppLoganEmis, SheppLoganBrainWeb)
        @test ImagePhantoms.shepp_logan_values(case()) isa Vector
        @test ellipse_parameters(case()) isa Matrix
    end
    @test ellipse_parameters(SouthPark()) isa Matrix

    image0 = @NOTinferred shepp_logan(2^6, SheppLoganEmis())
    @test image0 isa Matrix
    image1 = @NOTinferred shepp_logan(2^6, SheppLoganEmis(); oversample=1)
    image2 = @NOTinferred shepp_logan(2^6, SheppLoganEmis(); oversample=3)
    @test image0 == image2

    ob = shepp_logan(SheppLoganEmis())
    x = LinRange(-1,1,2^6+1) * 0.5
    y = LinRange(-1,1,2^6) * 0.5
    image = phantom(x, y, ob)
    @test image isa Matrix
    image = phantom(ob).(x,y')
    @test image isa Matrix

    r = LinRange(-0.5,0.5,2^5)
    ϕ = (0:30:180) * deg2rad(1)
    sino = radon(ob).(r,ϕ')
    @test sino isa Matrix
    sino = radon(r, ϕ, ob)
    @test sino isa Matrix

    kx = LinRange(-1,1,2^6) * 9
    ky = LinRange(-1,1,2^6+1) * 9
    kspace = spectrum(ob).(kx, ky')
    @test kspace isa Matrix
    kspace = spectrum(kx, ky, ob)
    @test kspace isa Matrix

    image = shepp_logan(40, 50, SouthPark(), fovs=(1,1))
end
