# mri-sense.jl test
# todo: move most to lit, simplify here

include("helper.jl") # todo

using MIRTjim: jim, prompt; jim(:prompt, true)
using ImagePhantoms: ellipse_parameters, SheppLoganBrainWeb, Ellipse
using ImagePhantoms: phantom, spectrum
using ImagePhantoms: smap_fit, spectra
using FFTW: fft, fftshift
using Unitful: mm
using Random: seed!
using MIRT: ir_mri_sensemap_sim
using ImageGeoms: embed
using LazyGrids: ndgrid
#using LinearAlgebra: norm
using Test: @test, @testset, @test_throws, @inferred

# define object
fovs = (256mm, 250mm)
nx,ny = (128, 100) .* 2
oa = ellipse_parameters(SheppLoganBrainWeb() ; disjoint=true, fovs)
seed!(0)
oa[:,end] = [1; 2; rand(ComplexF32,8)] # random phases
oa = Ellipse(oa)
	oa = oa[1:3] # todo simplify
dx, dy = fovs ./ (nx,ny)
x = (-(nx÷2):(nx÷2-1)) * dx
y = (-(ny÷2):(ny÷2-1)) * dy
oversample = 3
image1 = phantom(x, y, oa, oversample)
#jim(x, y, image1)

@test spectrum(oa[1]) isa Function

mask = trues(nx,ny)
mask[:,[1:2;end-2:end]] .= false
mask[[1:8;end-8:end],:] .= false
@assert mask .* image1 == image1
#jim(x, y, mask, "mask")


# smaps
ncoil = 2
#smap = ir_mri_sensemap_sim(dims=(nx,ny), ncoil=ncoil, orbit_start=[45])
#smap[:,:,1] .= 7 * exp.(2im * pi * (0:nx-1) / nx * 3)
smap = cat(dims=3,
	7 * exp.(im * 2pi * (0:nx-1) / nx * 3) * ones(ny)',
    5 * ones(nx) * exp.(-im * 2pi * (0:ny-1) / ny * 1)',
)
#=
=#
sfun = smap -> [real(smap);;;imag(smap)]
p1 = jim(x, y, sfun(smap), "Sensitivity maps raw"; prompt=false)

ssos = sqrt.(sum(abs.(smap).^2, dims=ndims(smap))) # SSoS
ssos = selectdim(ssos, ndims(smap), 1)
p2 = jim(ssos, "SSoS for ncoil=$ncoil"; prompt=false)

#= todo: no ssos normalization during tests!
for ic=1:ncoil
    selectdim(smap, ndims(smap), ic) ./= ssos
end
=#
smap .*= mask
p3 = jim(x, y, sfun(smap), "Sensitivity maps (masked)"; prompt=false)
#jim(p1, p2, p3)

#=
jim(
jim(real(smap), "real"),
jim(imag(smap)),
jim(abs.(smap), "abs"),
jim(angle.(smap)),
layout=(2,2),
)
=#

#smap = ones(ComplexF32, nx, ny)
#smap = cat(dims=3, [rand(ComplexF32)*smap for ic=1:ncoil]...)

#stacker = x -> [selectdim(x, ndims(x), i) for i=1:size(x)[end]]
stacker = x -> [(@view x[:,:,i]) for i=1:size(x,3)]
smaps = stacker(smap)
#jim(sfun(smaps), "smaps")


# fit each smap

deltas = (dx, dy)
kmax = 7

#=
plan = smap_basis(mask ; kmax, deltas)
coefs = [plan.B \ selectdim(smap, ndims(smap), ic)[mask] for ic=1:ncoil]
tmp = [embed(plan.B * coefs[ic], mask) for ic=1:ncoil]
tmp = cat(dims=3, tmp...)
#@show norm(tmp - smap) / norm(smap)
jim(tmp - smap)
=#

#=
function reala(v::Vector{<:AbstractArray})
	s = map(x -> ComplexF32.(x), v) # single precision
	r = map(real, s)
	if r ≈ s
		return r
	end
	i = map(imag, s)
	if i ≈ im * s
		return i
	end
	return s
end
=#

fit = smap_fit(smaps, embed; mask, kmax, deltas)
#smaps_fit = fit.smaps
#if map(real, tmp) ≈ tmp
#= for testing:
	fit = (B = fit.B, coefs = reala(fit.coefs), nrmse = fit.nrmse,
		smaps = reala(fit.smaps), ν = fit.ν)
=#
#fun = x -> real(ComplexF32.(x))
 #map(real, fit.smaps) ≈ fit.smaps
#map(x -> real(ComplexF32.(x)), fit.smaps) ≈ map(x -> ComplexF32.(x), fit.smaps
 #ComplexF32.(fit.smaps[1]) ≈ real(ComplexF32.(fit.smaps[1]))

coefs = map(x -> reshape(x,15,15), fit.coefs)
jim(
 jim(x, y, sfun(smaps), "orig", prompt=false),
 jim(x, y, sfun(fit.smaps), "fit", prompt=false),
 jim(x, y, sfun(1e6*(fit.smaps - smaps)), "err * 1e6", prompt=false),
 jim(-kmax:kmax, -kmax:kmax, sfun(coefs), "coefs", prompt=false),
#jim(-kmax:kmax, -kmax:kmax, real(coefs[1]), prompt=false),
#jim(-kmax:kmax, -kmax:kmax, imag(coefs[1]), prompt=false),
#jim(-kmax:kmax, -kmax:kmax, real(coefs[2]), prompt=false),
#jim(-kmax:kmax, -kmax:kmax, imag(coefs[2]), prompt=false),
)

fx = (-(nx÷2):(nx÷2-1)) / (nx*dx)
fy = (-(ny÷2):(ny÷2-1)) / (ny*dy)
gx, gy = ndgrid(fx, fy)


# test one object with no smap

ob = oa[3]
coil = 2

image1 = phantom(x, y, [ob], oversample)

fun0 = spectrum(ob)
#kspace0 = fun0.(vec(gx), vec(gy))
kspace0 = fun0.(fx, fy')
#kspace0 = reshape(kspace0, nx, ny)
kspace1 = spectrum(fx, fy, [ob])
@test kspace0 / oneunit(eltype(kspace0)) ≈ kspace1 / oneunit(eltype(kspace1))

kspace1 = myfft(image1) * dx * dy

p0 = jim(fx, fy, kspace0; prompt=false)
p1 = jim(fx, fy, kspace1; prompt=false)
p2 = jim(fx, fy, kspace1 - kspace0; prompt=false)
#jim(p0, p1, p2)


# test one object with one smap

image2 = image1 .* smaps[coil]
p1 = jim(x, y, sfun(image1); prompt=false)
p2 = jim(x, y, sfun(image2); prompt=false)
p3 = jim(x, y, sfun(smaps[coil]); prompt=false)
jim(p1, p2, p3)
throw(7)


kspace2 = myfft(image2) * dx * dy
fun = spectrum(ob, fit, coil)
#kspace3 = fun.(vec(gx), vec(gy))
kspace3 = -fun.(fx, fy') # todo: why negative!?
p0 = jim(fx, fy, kspace2, "fft"; prompt=false)
p1 = jim(fx, fy, kspace3, "anal"; prompt=false)
p2 = jim(fx, fy, kspace3 - kspace2, "err"; prompt=false)
jim(p0, p1, p2)

for i=1:100
	local fun = real
	jim(fx, fy, fun.(kspace1), "1")
	jim(fx, fy, fun.(kspace2), "2")
	jim(fx, fy, fun.(kspace3), "3")
end
#=
=#
throw(8)


#=
fun0 = spectrum(ob)
kspace0 = fun0.(vec(gx), vec(gy))
kspace0 = reshape(kspace0, nx, ny)
kspace1 = spectrum(fx, fy, [ob])
@test kspace0 / oneunit(eltype(kspace0)) ≈ kspace1 / oneunit(eltype(kspace1))




fun = spectrum(ob, fit, coil)
kspace3 = fun.(vec(gx), vec(gy))
kspace3 = reshape(kspace2, nx, ny)
p3 = jim(fx, fy, kspace1)
p4 = jim(fx, fy, kspace2)


#kspace1 = spectra(vec(gx), vec(gy), ob, fit)

# todo later: LinearMapAO for Array of Arrays?

#todo Xtrue = spectrum(-(nx÷2):(nx÷2-1), -(ny÷2):(ny÷2-1), object, 2)
#tmp = cat(dims=3, (fit.smaps - smaps)...) # todo: build into jim

=#
