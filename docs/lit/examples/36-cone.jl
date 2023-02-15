#=
# [Cone](@id 36-cone)

This page illustrates the `Cone` shape in the Julia package
[`ImagePhantoms`](https://github.com/JuliaImageRecon/ImagePhantoms.jl).
=#

#srcURL


# ### Setup

# Packages needed here.

using ImagePhantoms: Object, phantom, radon, spectrum
using ImagePhantoms: Cone, cone
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
is the cone,
specified by its center, base radii, height, angle(s) and value.
All of the methods in `ImagePhantoms` support physical units,
so we use such units throughout this example.
(Using units is recommended but not required.)

Here are 3 ways to define a `Object{Cone}`,
using physical units.
=#

center = (20mm, 10mm, -15mm)
width = (25mm, 30mm, 35mm) # x radius, y radius, height
ϕ0s = :(π/6) # symbol version for nice plot titles
ϕ0 = eval(ϕ0s)
angles = (ϕ0, 0, 0)
Object(Cone(), center, width, angles, 1.0f0) # top-level constructor
cone( 20mm, 10mm, 5mm, 25mm, 35mm, 15mm, π/6, 0, 0, 1.0f0) # 9 arguments
ob = cone(center, width, angles, 1.0f0) # tuples (recommended use)


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
   title="Cone, rotation ϕ=$ϕ0s", xlabel="x", ylabel="y")


# The image integral should match the object volume:
volume = IP.volume(ob)
(sum(img)*prod(ig.deltas), volume)


# Show middle slices
jim(mid3(img), "Middle 3 planes")
