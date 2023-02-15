#=
# [MRI SENSE](@id 10-mri-sense)

This page illustrates the `mri_smap_fit` and `mri_spectra` methods
in the Julia package
[`ImagePhantoms`](https://github.com/JuliaImageRecon/ImagePhantoms.jl)
for performing MRI simulations with realistic sensitivity encoding (SENSE).
=#

#srcURL

# ### Setup

# Packages needed here.

using ImagePhantoms: ellipse_parameters, SheppLoganBrainWeb, ellipse
using ImagePhantoms: phantom, mri_smap_fit, mri_spectra
using FFTW: fft, fftshift
using ImageGeoms: embed
using LazyGrids: ndgrid
using MIRTjim: jim, prompt
using Random: seed!
using Unitful: mm

# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);


#=
## Overview

Modern MRI scanners use multiple receive coils
each of which has its own "sensitivity map" (or "coil profile").
Realistic MRI simulations should account for the effects of those
sensitivity maps analytically, rather than committing the "inverse crime"
of using rasterized phantoms and maps.

See the 2012 paper
[Guerquin-Kern et al.](https://doi.org/10.1109/TMI.2011.2174158)
that combines analytical k-space values of the phantom
with an analytical model for the sensitivity maps.
This package follows the recommended approach from that paper.
We used the `mri_smap_fit` function to fit each sensitivity map
with a modest number of complex exponential signals.
Then, instead of using the `spectrum` function
we use the `spectra` function
to generate simulated k-space data
from analytical phantoms (like ellipses).
=#


# Because FFTW.fft cannot handle units, this function is a work-around.
function myfft(x::AbstractArray{T}) where {T <: Number}
    u = oneunit(T)
    return fftshift(fft(fftshift(x) / u)) * u
end


# ## Phantom

# Image geometry:

fovs = (256mm, 250mm)
nx, ny = (128, 100) .* 2
dx, dy = fovs ./ (nx,ny)
x = (-(nx÷2):(nx÷2-1)) * dx
y = (-(ny÷2):(ny÷2-1)) * dy

#=
Define Shepp-Logan phantom object,
with random complex phases
to make it a bit more realistic.
=#

params = ellipse_parameters(SheppLoganBrainWeb() ; disjoint=true, fovs)
seed!(0)
phases = [1; rand(ComplexF32,9)] # random phases
params = [(p[1:5]..., phases[i]) for (i, p) in enumerate(params)]
oa = ellipse(params)
oversample = 3
image0 = phantom(x, y, oa, oversample)
cfun = z -> cat(dims = ndims(z)+1, real(z), imag(z))
jim(x, y, cfun(image0), "Digital phantom\n (real | imag)")

#=
In practice, sensitivity maps are usually estimated
only over portion of the image array,
so we define a simple `mask` here
to exercise this issue.
=#

mask = trues(nx,ny)
mask[:,[1:2;end-2:end]] .= false
mask[[1:8;end-8:end],:] .= false
@assert mask .* image0 == image0
jim(x, y, mask, "mask")


#=
## Sensitivity maps

Here we use highly idealized sensitivity maps,
roughly corresponding to the
[Biot-Savart law](https://en.wikipedia.org/wiki/Biot-Savart_law)
for an infinite thin wire,
as a crude approximation of a
[birdcage coil](https://en.wikipedia.org/wiki/Radiofrequency_coil).
One wire is outside the upper right corner,
the other is outside the left border.
=#

# response at (x,y) to wire at (wx,wy)
function biot_savart_wire(x, y, wx, wy)
    phase = cis(atan(y-wy, x-wx))
    return oneunit(x) / sqrt(sum(abs2, (x-wx, y-wy))) * phase # 1/r falloff
end

ncoil = 2
wire1 = (a,b) -> biot_savart_wire(a, b, maximum(x) + 8dx, maximum(y) + 8dy)
wire2 = (a,b) -> biot_savart_wire(a, b, minimum(x) - 20dx, zero(dy))
smap = [wire1.(x, y'), wire2.(x, y')]
smap[1] *= cis(3π/4) # match coil phases at image center, ala "quadrature phase"
smap = cat(dims=3, smap...)
smap /= maximum(abs, smap)
mag = abs.(smap)
phase = angle.(smap)

jim(
 jim(x, y, mag, "|Sensitivity maps raw|"; color=:cividis, ncol=1, prompt=false),
 jim(x, y, phase, "∠(Sensitivity maps raw)"; color=:hsv, ncol=1, prompt=false),
)

#=
Typical sensitivity map estimation methods
normalize the maps
so that the square-root of the sum of squares (SSoS) is unity:
=#

ssos = sqrt.(sum(abs.(smap).^2, dims=ndims(smap))) # SSoS
ssos = selectdim(ssos, ndims(smap), 1)
jim(x, y, ssos, "SSoS for ncoil=$ncoil"; color=:cividis, clim=(0,1))

for ic=1:ncoil # normalize
    selectdim(smap, ndims(smap), ic) ./= ssos
end
smap .*= mask
stacker = x -> [(@view x[:,:,i]) for i=1:size(x,3)]
smaps = stacker(smap) # code hereafter expects vector of maps
jim(x, y, cfun(smaps), "Sensitivity maps (masked and normalized)")


#=
## Sensitivity map fitting using complex exponentials

The `mri_smap_fit` function fits each `smap`
with a linear combination of complex exponential signals.
(These signals are not orthogonal due to the `mask`.)
With frequencies `-9:9/N`, the maximum error is ≤ 0.4%.
=#

deltas = (dx, dy)
kmax = 9
fit = mri_smap_fit(smaps, embed; mask, kmax, deltas)
jim(
 jim(x, y, cfun(smaps), "Original maps"; prompt=false, clim=(-1,1)),
 jim(x, y, cfun(fit.smaps), "Fit maps"; prompt=false, clim=(-1,1)),
 jim(x, y, cfun(100 * (fit.smaps - smaps)), "error * 100"; prompt=false),
)

# The fit coefficients are smaller near `±kmax`
# so probably `kmax` is large enough.

coefs = map(x -> reshape(x, 2kmax+1, 2kmax+1), fit.coefs)
jim(-kmax:kmax, -kmax:kmax, cfun(coefs), "Coefficients")


# ## Compare FFT with analytical spectra

# Frequency sample vectors:
fx = (-(nx÷2):(nx÷2-1)) / (nx*dx) # crucial to match `mri_smap_basis` internals!
fy = (-(ny÷2):(ny÷2-1)) / (ny*dy)
gx, gy = ndgrid(fx, fy);

# Analytical spectra computation for complex phantom using all smaps.
# Note the `fit` argument.
kspace1 = mri_spectra(vec(gx), vec(gy), oa, fit)
kspace1 = [reshape(k, nx, ny) for k in kspace1]
p1 = jim(fx, fy, cfun(kspace1), "Analytical")

# FFT spectra computation based on digital image and sensitivity maps:
image2 = [image0 .* s for s in smaps] # digital
kspace2 = myfft.(image2) * (dx * dy)
p2 = jim(fx, fy, cfun(kspace2), "FFT-based")

#
p3 = jim(fx, fy, cfun(kspace2 - kspace1), "Error")


# Zoom in to illustrate similarity:

xlims = (-1,1) .* (0.06/mm)
ylims = (-1,1) .* (0.06/mm)
jim(
 jim(fx, fy, real(kspace1[1]), "Analytical"; xlims, ylims, prompt=false),
 jim(fx, fy, real(kspace2[1]), "FFT-based"; xlims, ylims, prompt=false),
 jim(fx, fy, real(kspace2[1] - kspace1[1]), "Error"; xlims, ylims, prompt=false),
)


#=
In summary, the `mri_smap_fit` and `mri_spectra` methods here
reproduce the approach in the 2012 Guerquin-Kern paper, cited above,
enabling parallel MRI simulations that avoid an inverse crime.
=#
