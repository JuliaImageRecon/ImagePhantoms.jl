# test/iter.jl

using ImagePhantoms: radon, rect, cuboid
using ImagePhantoms: _radon
using Unitful: m
using Test: @test, @testset, @inferred


@testset "iter2" begin
    ob = rect((4m, 3m), (2m, 5m), π/6, 1.0f0 / m^2)
    nr, dr = 2^4, 0.02m
    r = (-nr÷2:nr÷2-1) * dr
    ϕ = 0:3
    sino1 = @inferred radon(r, ϕ, [ob])
    itr = Iterators.product(r, ϕ)
    sino2 = @inferred radon(itr, [ob])
    @test sino1 == sino2
end

@testset "iter3" begin
    ob = cuboid((1m, 2m, 3m), (4m, 5m, 6m), (π/6, 0, 0), 1.0f0 / m^2)
    nu, du = 2^4, 0.02m
    nv, dv = 2^3, 0.03m
    u = (-nu÷2:nu÷2-1) * du
    v = (-nv÷2:nv÷2-1) * dv
    ϕ = 0:3
    θ = deg2rad.([0, 10])
    proj1 = @inferred radon(u, v, ϕ, θ, [ob])
    itr = Iterators.product(u, v, ϕ, θ)
    proj2 = @inferred radon(itr, [ob])
    @test proj1 == proj2
end
