#---------------------------------------------------------
# # [Ellipse](@id 2-ellipse)
#---------------------------------------------------------

# This page illustrates the `Ellipse` shape in the Julia package
# [`ImagePhantoms`](https://github.com/JuliaImageRecon/ImagePhantoms.jl).

# ### Setup

# Packages needed here.

using ImagePhantoms: Ellipse, phantom, radon, spectrum
using ImageGeoms: ImageGeom, axesf
using MIRTjim: jim, prompt
using FFTW: fft, fftshift
using Unitful: mm, unit
using UnitfulRecipes
using Plots: plot, plot!, scatter!, default; default(markerstrokecolor=:auto)

# The following line is helpful when running this example.jl file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);

# ### Overview

# One of the most basic shapes used in constructing 2D digital image phantoms
# is the ellipse, specified by its center, radii, angle and value.
# All of the methods in `ImagePhantoms` support physical units,
# so we use such units throughout this example.
# (Using units is recommended but not required.)

# Define an ellipse object, using physical units.
width = (2mm, 8mm)
ob = Ellipse((4mm, 3mm), width, π/6, 1.0f0)


# ### Phantom image using `phantom`

# Make a digital image of it using `phantom` and display it.
dx = 0.02mm
dy = 0.025mm
(M,N) = (2^10,2^10+2)
x = (-M÷2:M÷2-1) * dx
y = (-N÷2:N÷2-1) * dy
img = phantom(x, y, [ob])
jim(x, y, img)


# Hereafter we use `ImageGeoms` to simplify the indexing.

ig = ImageGeom(dims=(M,N), deltas=(dx,dy), offsets=(0.5,0.5))
@assert all(axes(ig) .≈ (x,y))
p1 = jim(axes(ig)..., img, "Ellipse phantom", xlabel="x", ylabel="y")


# ### Spectrum using `spectrum`

# Let's examine the spectrum of this image.
# There are two ways to do this:
# * using the analytical Fourier transform of the ellipse via `spectrum`
# * applying the DFT via FFT to the digital image.
# Because the shape has units `mm`, the spectra axes have units cycles/mm.

zscale = 4 / π / prod(width) # normalize spectra by area
spectrum_exact = spectrum(axesf(ig)..., [ob]) * zscale
sp = z -> max(log10(abs(z)/oneunit(abs(z))), -6) # log-scale for display
clim = (-6, 0) # colorbar limit for display
(xlabel, ylabel) = ("ν₁", "ν₂")
p2 = jim(axesf(ig)..., sp.(spectrum_exact), "log10|Spectrum|"; clim, xlabel, ylabel)


# Sadly `fft` cannot handle units currently, so this function is a work-around:
function myfft(x)
    u = unit(eltype(x))
    return fftshift(fft(fftshift(x) / u)) * u
end

# fx = (-M÷2:M÷2-1) / M / dx # appropriate frequency axes for DFT,
# fy = (-N÷2:N÷2-1) / N / dy # that are provided by axesf(ig)
spectrum_fft = myfft(img) * dx * dy * zscale
p3 = jim(axesf(ig)..., sp.(spectrum_fft), "log10|DFT|"; clim, xlabel, ylabel)


# Compare the DFT and analytical spectra to validate the code
@assert maximum(abs, spectrum_exact - spectrum_fft) /
        maximum(abs, spectrum_exact) < 2e-2
p4 = jim(axesf(ig)..., abs.(spectrum_fft - spectrum_exact), "Difference"; xlabel, ylabel)
jim(p1, p4, p2, p3)


# ### Radon transform using `radon`

# Examine the Radon transform of the ellipse using `radon`,
# and validate it using the projection-slice theorem aka Fourier-slice theorem.

dr = 0.02mm # radial sample spacing
nr = 2^10 # radial sinogram bins
r = (-nr÷2:nr÷2-1) * dr # radial samples
fr = (-nr÷2:nr÷2-1) / nr / dr # corresponding spectral axis
ϕ = deg2rad.(0:180) # * Unitful.rad # todo round unitful Unitful.°
sino = radon(ob).(r, ϕ') # sample Radon transform of a single shape object
p5 = jim(r, rad2deg.(ϕ), sino; aspect_ratio=:none, title="sinogram", yflip=false, xlabel="r", ylabel="ϕ")

#src clim=(0, 2*maximum(width)*ob.value), # todo

# Note that the maximum sinogram value is about 16mm which makes sense
# for an ellipse whose long axis has "radius" 8mm.

# The above sampling generated a parallel-beam sinogram,
# but one could make a fan-beam sinogram simply by sampling `(r, ϕ)` appropriately.


# ### Fourier-slice theorem illustration

# Pick one particular view angle (55°) and look at its slice and spectra.
ia = argmin(abs.(ϕ .- deg2rad(55)))
slice = sino[:,ia]
slice_fft = myfft(slice) * dr
angle = round(rad2deg(ϕ[ia]), digits=1)

kx, ky = (fr * cos(ϕ[ia]), fr * sin(ϕ[ia])) # Fourier-slice theorem
slice_ft = spectrum(ob).(kx, ky)
@assert maximum(abs, slice_ft - slice_fft) / maximum(abs, slice_ft) < 2e-4

p3 = plot(r, slice, title="profile at ϕ = $angle", label="")
p4 = plot(title="1D spectra")
scatter!(fr, abs.(slice_fft), label="abs fft", color=:blue)
scatter!(fr, real(slice_fft), label="real fft", color=:green)
scatter!(fr, imag(slice_fft), label="imag fft", color=:red, xlims=(-1,1).*(1.0/mm))

plot!(fr, abs.(slice_ft), label="abs", color=:blue)
plot!(fr, real(slice_ft), label="real", color=:green)
plot!(fr, imag(slice_ft), label="imag", color=:red)
plot(p1, p5, p3, p4)


# The good agreement between the analytical spectra (solid lines)
# and the DFT samples (disks)
# validates that `phantom`, `radon`, and `spectrum`
# are all self consistent for this `Ellipse` object.
