#=
# [Cylinder](@id 35-cylinder)

This page illustrates the `Cylinder` shape in the Julia package
[`ImagePhantoms`](https://github.com/JuliaImageRecon/ImagePhantoms.jl).
=#

#srcURL


# ### Setup

# Packages needed here.

using ImagePhantoms: Object, phantom, radon, spectrum
using ImagePhantoms: Cylinder, cylinder
import ImagePhantoms as IP
using ImageGeoms: ImageGeom, axesf
using MIRTjim: jim, prompt, mid3
using FFTW: fft, fftshift, ifftshift
using LazyGrids: ndgrid
using Unitful: mm, unit, °
using Plots: plot, plot!, scatter!, default
using Plots # gif @animate
default(markerstrokecolor=:auto)


# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);


#=
## Overview

A basic shape used in constructing 3D digital image phantoms
is the cylinder,
specified by its center, radii, height, angle(s) and value.
All of the methods in `ImagePhantoms` support physical units,
so we use such units throughout this example.
(Using units is recommended but not required.)

Here are 3 ways to define a `Object{Cylinder}`,
using physical units.
=#

center = (20mm, 10mm, 5mm)
width = (25mm, 35mm, 15mm) # x radius, y radius, height
ϕ0s = :(π/6) # symbol version for nice plot titles
ϕ0 = eval(ϕ0s)
angles = (ϕ0, 0, 0)
Object(Cylinder(), center, width, angles, 1.0f0) # top-level constructor
cylinder( 20mm, 10mm, 5mm, 25mm, 35mm, 15mm, π/6, 0, 0, 1.0f0) # 9 arguments
ob = cylinder(center, width, angles, 1.0f0) # tuples (recommended use)


#=
## Phantom image using `phantom`

Make a 3D digital image of it using `phantom` and display it.
We use `ImageGeoms` to simplify the indexing.
=#

deltas = (1.0mm, 1.1mm, 0.9mm)
dims = (2^8, 2^8+2, 49) # odd
ig = ImageGeom( ; dims, deltas, offsets=:dsp)
oversample = 3
img = phantom(axes(ig)..., [ob], oversample)
p1 = jim(axes(ig), img;
   title="Cylinder, rotation ϕ=$ϕ0s", xlabel="x", ylabel="y")


# The image integral should match the object volume:
volume = IP.volume(ob)
(sum(img)*prod(ig.deltas), volume)


# Show middle slices
jim(mid3(img), "Middle 3 planes")


#=
## Spectrum using `spectrum`

There are two ways to examine the spectrum of this 3D image:
* using the analytical Fourier transform of the object via `spectrum`
* applying the DFT via FFT to the digital image.
Because the shape has units `mm`, the spectra axes have units cycles/mm.
Appropriate frequency axes for DFT are provided by `axesf(ig)`.
=#

vscale = 1 / volume # normalize spectra by volume
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
@assert errf < 4e-2
p4 = jim(axesf(ig), 1e3*abs.(spectrum_fft - spectrum_exact);
   title="|Difference| × 10³", xlabel, ylabel)
jim(p1, p4, p2, p3)


#=
## Parallel-beam projections using `radon`

Compute 2D projection views of the object using `radon`.
Validate it using the projection-slice theorem aka Fourier-slice theorem.
=#

pg = ImageGeom((2^8,2^7), (0.6mm,1.0mm), (0.5,0.5)) # projection sampling
ϕs, θs = (:(π/3), ϕ0s), (:(π/7), :(0))
ϕ3 = ϕ0
θ3 = atan(ob.width[3]/2, maximum(ob.width[1:2]))
ϕ, θ = [eval.(ϕs)..., ϕ3], [eval.(θs)..., θ3]
proj3 = [radon(axes(pg)..., ϕ[i], θ[i], [ob]) for i in 1:3] # 3 projections
smax = ob.value * sqrt(maximum(abs2, 2 .* ob.width[1:2]) + abs2(ob.width[3]))
p5 = jim(axes(pg)..., proj3; xlabel="u", ylabel="v", nrow = 1, title =
    "Projections at (ϕ,θ) = ($(ϕs[1]), $(θs[1])) and ($(ϕs[2]), $(θs[2]))\n
    and along long axis")


#=
Because the cylinder has maximum diameter of 70mm and height 15mm,
and one of the two views above was along the corresponding axis,
the maximum projection value is about
`sqrt(70^2 + 15^2)` ≈ 71.6mm.
=#

maxes = round.((smax, maximum.(proj3)...) ./ 1mm; digits=2)


#=
The integral of each projection should match the object volume:
=#

vols = round.(((p -> sum(p)*prod(pg.deltas)).(proj3)..., volume) ./ 1mm^3; digits=2)


# Look at a set of projections as the views orbit around the object.
ϕd = 0:6:360
ϕs = deg2rad.(ϕd)
θs = :(π/7)
θ = eval(θs)
projs = radon(axes(pg)..., ϕs, [θ], [ob]) # many projection views

if isinteractive()
    jim(axes(pg)..., projs; title="projection views $(ϕd)")
else
    anim = @animate for ip in 1:length(ϕd)
        jim(axes(pg), projs[:,:,ip,1]; xlabel="u", ylabel="v", prompt=false,
            title="ϕ=$(ϕd[ip])° θ=$θs", clim = (0,1) .* smax)
    end
    gif(anim, "cylinder.gif", fps = 6)
end


#=
The above sampling generated a parallel-beam projection,
but one could make a cone-beam projection
by sampling `(u, v, ϕ, θ)` appropriately.
See Sinograms.jl.
=#


#=
## Fourier-slice theorem illustration

Pick one particular view and compare its FFT
to a slice through the 3D object spectrum.
=#

ϕs, θs = :(π/3), :(π/7)
ϕ, θ = eval.((ϕs, θs))
proj = radon(axes(pg)..., ϕ, θ, [ob])
p6 = jim(axes(pg), proj; xlabel="u", ylabel="v", prompt=false,
    title = "Projection at (ϕ,θ) = ($ϕs, $θs)")

e1 = (cos(ϕ), sin(ϕ), 0)
e3 = (sin(ϕ)*sin(θ), -cos(ϕ)*sin(θ), cos(θ))
fu, fv = ndgrid(axesf(pg)...)
ff = vec(fu) * [e1...]' + vec(fv) * [e3...]' # fx,fy,fz for Fourier-slice theorem
spectrum_slice = spectrum(ob).(ff[:,1], ff[:,2], ff[:,3]) * vscale
spectrum_slice = reshape(spectrum_slice, pg.dims)
clim = (-6, 0) # colorbar limit for display
(xlabel, ylabel) = ("νᵤ", "νᵥ")
p7 = jim(axesf(pg), sp.(spectrum_slice); prompt=false,
    title = "log10|Spectrum Slice|", clim, xlabel, ylabel)
proj_fft = myfft(proj) * prod(pg.deltas) * vscale
p8 = jim(axesf(pg), sp.(proj_fft); prompt=false,
     title = "log10|FFT Spectrum|", clim, xlabel, ylabel)

errs = maximum(abs, spectrum_slice - proj_fft) / maximum(abs, spectrum_slice)
@assert errs < 1e-3
p9 = jim(axesf(pg), 1e3*abs.(proj_fft - spectrum_slice);
    title="|Difference| × 10³", xlabel, ylabel, prompt=false)
jim(p6, p7, p8, p9)


#=
The good agreement between
the 2D slice through the 3D analytical spectrum
and the FFT of the 2D projection view
validates that `phantom`, `radon`, and `spectrum`
are all self consistent
for this shape.
=#
