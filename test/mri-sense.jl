# mri-sense.jl test

using ImagePhantoms: ellipse_parameters, SheppLoganBrainWeb, ellipse
using ImagePhantoms: phantom, spectrum
using ImagePhantoms: mri_smap_fit, mri_spectra
using FFTW: fft, fftshift
using Unitful: mm
using Random: seed!
using ImageGeoms: embed
using LazyGrids: ndgrid
using LinearAlgebra: norm
using Test: @test, @testset, @test_throws, @inferred

# image geometry
fovs = (256mm, 250mm)
nx, ny = (64, 50) .* 1
dx, dy = fovs ./ (nx,ny)
x = (-(nx÷2):(nx÷2-1)) * dx
y = (-(ny÷2):(ny÷2-1)) * dy

# define object
oa = ellipse_parameters(SheppLoganBrainWeb() ; disjoint=true, fovs)
seed!(0)
oa[:,end] = [1; rand(ComplexF32,9)] # random phases
oa = ellipse(oa)
oversample = 3
image0 = phantom(x, y, oa, oversample)
cfun = z -> cat(dims = ndims(z)+1, real(z), imag(z))

@test spectrum(oa[1]) isa Function

# define mask
mask = trues(nx,ny)
mask[:,1:2] .= false
mask[[1:8;end-8:end],:] .= false
@assert mask .* image0 == image0

# sensitivity maps (smaps) (not normalized)
ncoil = 2
smap = cat(dims=3,
    7 * exp.(im * 2pi * ((0:nx-1) .- 27) / nx * 3) * ones(ny)',
    5 * ones(nx) * exp.(-im * 2pi * ((0:ny-1) .+ 39) / ny * 1)',
)
smap .*= mask
stacker = x -> [(@view x[:,:,i]) for i=1:size(x,3)]
smaps = stacker(smap)


# fit each smap
deltas = (dx, dy)
kmax = 7
fit = @NOTinferred mri_smap_fit(smaps, embed; mask, kmax, deltas)
@test fit isa NamedTuple
@test fit.nrmse ≤ 1e-6

coefs = map(x -> reshape(x,15,15), fit.coefs)

fx = (-(nx÷2):(nx÷2-1)) / (nx*dx) # crucial to match smap_basis internals!
fy = (-(ny÷2):(ny÷2-1)) / (ny*dy)
gx, gy = ndgrid(fx, fy)


# test one object with no smap
ob = oa[3]
coil = 2
image1 = phantom(x, y, [ob], oversample)

#=

fun0 = @inferred spectrum(ob)
kspace0 = fun0.(fx, fy')
kspace1 = spectrum(fx, fy, [ob])
#@test kspace0 ≈ kspace1 # https://github.com/PainterQubits/Unitful.jl/pull/468
@test kspace0 / oneunit(eltype(kspace0)) ≈ kspace1 / oneunit(eltype(kspace1))

kspace1 = myfft(image1) * dx * dy

#@test norm
=#


# test one object with one smap

image2 = image1 .* smaps[coil] # digital

kspace2 = myfft(image2) * dx * dy
fun = spectrum(ob, fit, coil)
kspace3 = fun.(fx, fy')

@test norm(kspace3 - kspace2) / norm(kspace2) ≤ 0.09


# test multiple objects with all smaps

#kspace0 = spectrum(oa).(vec(gx), vec(gy))
#kspace0 = reshape(spectrum(oa).(vec(gx), vec(gy)), nx, ny)
#kspace0 = spectrum(oa).(fx, fy') # without smaps
#kspace1 = mri_spectra(fx, fy, oa, fit) # no!
kspace1 = mri_spectra(vec(gx), vec(gy), oa, fit)
kspace1 = [reshape(k, nx, ny) for k in kspace1]
#@test kspace0 / oneunit(eltype(kspace0)) ≈ kspace1 / oneunit(eltype(kspace1[1]))

image2 = [image0 .* s for s in smaps] # digital
kspace2 = myfft.(image2) * (dx * dy)

norma = za -> sqrt(sum(abs2, norm.(za))) # for array of arrays
@test norma(kspace2 - kspace1) / norma(kspace2) ≤ 0.03
