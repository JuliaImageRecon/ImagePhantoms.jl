#---------------------------------------------------------
# # [Ellipse](@id 02-ellipse)
#---------------------------------------------------------

#=
This page illustrates the `Ellipse` shape in the Julia package
[`ImagePhantoms`](https://github.com/JuliaImageRecon/ImagePhantoms.jl).

This page was generated from a single Julia file:
[02-ellipse.jl](@__REPO_ROOT_URL__/02-ellipse.jl).
=#

#md # In any such Julia documentation,
#md # you can access the source code
#md # using the "Edit on GitHub" link in the top right.

#md # The corresponding notebook can be viewed in
#md # [nbviewer](http://nbviewer.jupyter.org/) here:
#md # [`02-ellipse.ipynb`](@__NBVIEWER_ROOT_URL__/02-ellipse.ipynb),
#md # and opened in [binder](https://mybinder.org/) here:
#md # [`02-ellipse.ipynb`](@__BINDER_ROOT_URL__/02-ellipse.ipynb).


# ### Setup

# Packages needed here.

using ImagePhantoms: ellipse
using ImagePhantoms: phantom, radon, spectrum
import ImagePhantoms as IP
using ImageGeoms: ImageGeom, axesf
using MIRTjim: jim, prompt
using FFTW: fft, fftshift, ifftshift
using Unitful: mm, unit, °
using Plots: plot, plot!, scatter!, default
default(markerstrokecolor=:auto)

# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);


# ### Overview

#=
A useful shape
for constructing 2D digital image phantoms
is the ellipse,
specified by its center, radii, angle and value.
All of the methods in `ImagePhantoms` support physical units,
so we use such units throughout this example.
(Using units is recommended but not required.)

Define an ellipse object,
using physical units.
=#

radii = (20mm, 80mm)
ob = ellipse((40mm, 30mm), radii, π/6, 1.0f0)


# ### Phantom image using `phantom`

# Make a digital image of it using `phantom` and display it.
dx, dy = 0.8mm, 1.0mm
M, N = (2^8, 2^8+2)
x = (-M÷2:M÷2-1) * dx
y = (-N÷2:N÷2-1) * dy
oversample = 2
img = phantom(x, y, [ob], oversample)
jim(x, y, img, "Ellipse image")


# Hereafter we use `ImageGeoms` to simplify the indexing.

M, N = (2^8, 2^8+17) # odd
ig = ImageGeom(dims=(M,N), deltas=(dx,dy), offsets=:dsp)
@assert axes(ig)[1] ≈ x
oversample = 2
img = phantom(axes(ig)..., [ob], oversample)
p1 = jim(axes(ig), img, "Ellipse phantom", xlabel="x", ylabel="y")


# The image integral should approximate the object area
area = IP.area(ob)
(sum(img) * prod(ig.deltas), area)


#=
### Spectrum using `spectrum`

There are two ways to examine the spectrum of this image:
* using the analytical Fourier transform of the object via `spectrum`
* applying the DFT via FFT to the digital image.
Because the shape has units `mm`, the spectra axes have units cycles/mm.
=#

vscale = 1 / area # normalize spectra by area
spectrum_exact = spectrum(axesf(ig)..., [ob]) * vscale
sp = z -> max(log10(abs(z)/oneunit(abs(z))), -6) # log-scale for display
clim = (-6, 0) # colorbar limit for display
(xlabel, ylabel) = ("ν₁", "ν₂")
p2 = jim(axesf(ig), sp.(spectrum_exact), "log10|Spectrum|"; clim, xlabel, ylabel)


# Sadly `fft` cannot handle units currently, so this function is a work-around:
function myfft(x)
    u = unit(eltype(x))
    return fftshift(fft(ifftshift(x) / u)) * u
end

spectrum_fft = myfft(img) * (prod(ig.deltas) * vscale)
p3 = jim(axesf(ig), sp.(spectrum_fft), "log10|DFT|"; clim, xlabel, ylabel)


# Compare the DFT and analytical spectra to validate the code
errf = maximum(abs, spectrum_exact - spectrum_fft) / maximum(abs, spectrum_exact)
@assert errf < 2e-3
p4 = jim(axesf(ig), 1e3*abs.(spectrum_fft - spectrum_exact),
    "|Difference| × 10³"; xlabel, ylabel)
jim(p1, p4, p2, p3)


# ### Radon transform using `radon`

# Examine the Radon transform of the object using `radon`,
# and validate it using the projection-slice theorem aka Fourier-slice theorem.

dr = 0.2mm # radial sample spacing
nr = 2^10 # radial sinogram bins
r = (-nr÷2:nr÷2-1) * dr # radial samples
fr = (-nr÷2:nr÷2-1) / nr / dr # corresponding spectral axis
ϕ = deg2rad.(0:180)
sino = radon(ob).(r, ϕ') # sample Radon transform of a single shape object
smax = ob.value * 2 * maximum(radii)
p5 = jim(r, rad2deg.(ϕ), sino; title="sinogram",
    xlabel="r", ylabel="ϕ", clim = (0,1) .* smax)


#=
The maximum sinogram value is about
160mm,
which makes sense for
an ellipse whose long axis has "radius" 80mm.

The above sampling generated a parallel-beam sinogram,
but one could make a fan-beam sinogram by sampling `(r, ϕ)` appropriately.
=#


# ### Fourier-slice theorem illustration

# Pick one particular view angle (55°) and look at its slice and spectra.
ia = argmin(abs.(ϕ .- 55°))
slice = sino[:,ia]
slice_fft = myfft(slice) * dr
ϕd = round(rad2deg(ϕ[ia]), digits=1)

kx, ky = (fr * cos(ϕ[ia]), fr * sin(ϕ[ia])) # Fourier-slice theorem
slice_ft = spectrum(ob).(kx, ky)
errs = maximum(abs, slice_ft - slice_fft) / maximum(abs, slice_ft)
@assert errs < 2e-4

p3 = plot(r, slice, title="profile at ϕ = $ϕd", label="")
p4 = plot(title="1D spectra")
scatter!(fr, abs.(slice_fft), label="abs fft", color=:blue)
scatter!(fr, real(slice_fft), label="real fft", color=:green)
scatter!(fr, imag(slice_fft), label="imag fft", color=:red,
    xlims=(-1,1) .* (0.1/mm))

plot!(fr, abs.(slice_ft), label="abs", color=:blue)
plot!(fr, real(slice_ft), label="real", color=:green)
plot!(fr, imag(slice_ft), label="imag", color=:red)
plot(p1, p5, p3, p4)


#=
The good agreement between the analytical spectra (solid lines)
and the DFT samples (disks)
validates that `phantom`, `radon`, and `spectrum`
are all self consistent
for this shape.
=#
