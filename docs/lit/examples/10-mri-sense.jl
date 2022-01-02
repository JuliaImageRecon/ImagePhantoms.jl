#---------------------------------------------------------
# # [MRI SENSE](@id 10-mri-sense)
#---------------------------------------------------------

#=
This page illustrates the `mri_smap_fit` and `mri_spectra` methods
in the Julia package
[`ImagePhantoms`](https://github.com/JuliaImageRecon/ImagePhantoms.jl)
for performing MRI simulations with realistic sensitivity encoding (SENSE).

This page was generated from a single Julia file:
[10-mri-sense.jl](@__REPO_ROOT_URL__/10-mri-sense.jl).
=#
#md # In any such Julia documentation,
#md # you can access the source code
#md # using the "Edit on GitHub" link in the top right.

#md # The corresponding notebook can be viewed in
#md # [nbviewer](http://nbviewer.jupyter.org/) here:
#md # [`10-mri-sense.ipynb`](@__NBVIEWER_ROOT_URL__/mri/10-mri-sense.ipynb),
#md # and opened in [binder](https://mybinder.org/) here:
#md # [`10-mri-sense.ipynb`](@__BINDER_ROOT_URL__/mri/10-mri-sense.ipynb).


# ### Setup

# Packages needed here.

using ImagePhantoms: ellipse_parameters, SheppLoganBrainWeb, Ellipse
using ImagePhantoms: phantom, mri_smap_fit, mri_spectra
using FFTW: fft, fftshift
using ImageGeoms: embed
using LazyGrids: ndgrid
using MIRT: ir_mri_sensemap_sim
using MIRTjim: jim, prompt; jim(:prompt, true)
using Random: seed!
using Unitful: mm

# The following line is helpful when running this example.jl file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);


# ### Overview

#=
Modern MRI scanners use multiple receive coils
each of which has its own "sensitivity map" (or "coil profile").
Realistic MRI simulations should account for the effects of those
sensitivity maps analytically, rather than committing the "inverse crime"
of using rasterized phantoms and maps.

See the 2012 paper
[Guerquin-Kern et al.](http://doi.org/10.1109/TMI.2011.2174158)
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


# ### Phantom

# Image geometry:

fovs = (256mm, 250mm)
nx, ny = (128, 100) .* 2
dx, dy = fovs ./ (nx,ny)
x = (-(nx÷2):(nx÷2-1)) * dx
y = (-(ny÷2):(ny÷2-1)) * dy

# Define Shepp-Logan phantom object, with random complex phases
# to make it a bit more realistic.

oa = ellipse_parameters(SheppLoganBrainWeb() ; disjoint=true, fovs)
seed!(0)
oa[:,end] = [1; randn(ComplexF32, 9)] # random phases
oa = Ellipse(oa)
oversample = 3
image0 = phantom(x, y, oa, oversample)
cfun = z -> cat(dims = ndims(z)+1, real(z), imag(z))
jim(x, y, cfun(image0), "Digital phantom\n (real | imag)"; ncol=1)

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


# ### Sensitivity maps

# Here we use simulated sensitivity maps from MIRT:

ncoil = 2
smap = ir_mri_sensemap_sim(dims=(nx,ny), ncoil=ncoil, orbit_start=[0])
jim(x, y, cfun(smap), "Sensitivity maps raw")

#=
Typical sensitivity map estimation methods
normalize the maps
so that the square-root of the sum of squares (SSoS) is unity:
=#

ssos = sqrt.(sum(abs.(smap).^2, dims=ndims(smap))) # SSoS
ssos = selectdim(ssos, ndims(smap), 1)
jim(ssos, "SSoS for ncoil=$ncoil")

for ic=1:ncoil # normalize
    selectdim(smap, ndims(smap), ic) ./= ssos
end
smap .*= mask
stacker = x -> [(@view x[:,:,i]) for i=1:size(x,3)]
smaps = stacker(smap) # code hereafter expects vector of maps
jim(x, y, cfun(smaps), "Sensitivity maps (masked and normalized)")


# ### Sensitivity map fitting using complex exponentials

#=
The `mri_smap_fit` function fits each `smap`
with a linear combination of complex exponential signals.
(These signals are not orthogonal due to the `mask`.)
With frequencies `-7:7/N`, the maximum error is ≤ 0.2%.
=#

deltas = (dx, dy)
kmax = 7
fit = mri_smap_fit(smaps, embed; mask, kmax, deltas)
jim(
 jim(x, y, cfun(smaps), "Original maps"; prompt=false, clim=(-1,1)),
 jim(x, y, cfun(fit.smaps), "Fit maps"; prompt=false, clim=(-1,1)),
 jim(x, y, cfun(100 * (fit.smaps - smaps)), "error * 100"; prompt=false),
)

# The fit coefficients are smaller near `kmax`
# so probably `kmax` is large enough.

coefs = map(x -> reshape(x,15,15), fit.coefs)
jim(-kmax:kmax, -kmax:kmax, cfun(coefs), "coefs", prompt=false)


# ### Compare FFT with analytical spectra

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
 jim(fx, fy, real(kspace1[1]), "Analytical"; xlims, ylims),
 jim(fx, fy, real(kspace2[1]), "FFT-based"; xlims, ylims),
 jim(fx, fy, real(kspace2[1] - kspace1[1]), "Error"; xlims, ylims),
)


#=
In summary, the `mri_smap_fit` and `mri_spectra` methods here
reproduce the approach in the 2012 Guerquin-Kern paper, cited above,
enabling parallel MRI simulations that avoid an inverse crime.
=#
