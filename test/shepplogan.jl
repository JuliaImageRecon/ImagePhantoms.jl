#=
shepplogan.jl
=#

#using MIRTjim: jim, prompt
using ImagePhantoms # many
using Test: @test, @testset, @test_throws, @inferred


@testset "shepp" begin
    for case in (SheppLoganToft, SheppLoganEmis, SheppLoganBrainWeb, SouthPark)
        @test ellipse_parameters(case()) isa Matrix
    end

    image = shepp_logan(256, SheppLoganEmis())
    @test image isa Matrix
#   jim(image)

    ob = shepp_logan(SheppLoganEmis())
    x = LinRange(-1,1,201) * 0.5
    y = LinRange(-1,1,200) * 0.5
    image = phantom(x, y, ob)
    @test image isa Matrix
    image = phantom(ob).(x,y')
    @test image isa Matrix
#   jim(x, y, image)

    r = LinRange(-0.5,0.5,101)
    ϕ = deg2rad.(0:180)
    sino = radon(ob).(r,ϕ')
    @test sino isa Matrix
    sino = radon(r, ϕ, ob)
    @test sino isa Matrix
#   jim(r, ϕ, sino; aspect_ratio=:none, yflip=false)

    kx = LinRange(-1,1,100) * 9
    ky = LinRange(-1,1,101) * 9
    kspace = spectrum(ob).(kx, ky')
    @test kspace isa Matrix
    kspace = spectrum(kx, ky, ob)
    @test kspace isa Matrix
#   jim(kx, ky, kspace)

    image = shepp_logan(80, 100, SouthPark(), fovs=(1,1))
#   jim(image; yflip=false)
end
