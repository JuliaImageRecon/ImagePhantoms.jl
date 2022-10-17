#---------------------------------------------------------
# # [ImagePhantoms overview](@id 01-overview)
#---------------------------------------------------------

#=
This page explains the Julia package
[`ImagePhantoms`](https://github.com/JuliaImageRecon/ImagePhantoms.jl).

This page was generated from a single Julia file:
[01-overview.jl](@__REPO_ROOT_URL__/01-overview.jl).
=#

#md # In any such Julia documentation,
#md # you can access the source code
#md # using the "Edit on GitHub" link in the top right.

#md # The corresponding notebook can be viewed in
#md # [nbviewer](http://nbviewer.jupyter.org/) here:
#md # [`01-overview.ipynb`](@__NBVIEWER_ROOT_URL__/01-overview.ipynb),
#md # and opened in [binder](https://mybinder.org/) here:
#md # [`01-overview.ipynb`](@__BINDER_ROOT_URL__/01-overview.ipynb).


# ### Setup

# Packages needed here.

#nb import Pkg
#nb Pkg.add(["ImageGeoms", "MIRTjim", "Unitful"])

using ImagePhantoms
using ImageGeoms: ImageGeom, axesf
using MIRTjim: jim, prompt
using Plots; default(markerstrokecolor=:auto, label="")
using Unitful: mm
using InteractiveUtils: versioninfo

# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);


# ### Overview

#=
When developing image reconstruction methods,
it can be helpful to simulate data (e.g., sinograms)
using software-defined images called phantoms.

The simplest method here is to make a Shepp-Logan phantom image
similar its use in other packages.
=#

image = shepp_logan(256) # CT version by default
jim(image, "SheppLogan"; clim=(0.9, 1.1))


# ## Sinograms

#=
Often for image reconstruction algorithm development,
we need not only the phantom image, but also its
[sinogram, i.e., Radon transform](https://en.wikipedia.org/wiki/Radon_transform)
and spectrum.

The 2D sinogram coordinate system used here is
```math
p(r, ϕ) = \int_{-∞}^{∞} f(r \cos ϕ - ℓ \sin ϕ, r \sin ϕ + ℓ \cos ϕ) \, \mathrm{d} ℓ.
```

We start with the vector of ellipses that defines the phantom,
using a typical field of view (FOV) of 200mm for a head:
=#

objects = shepp_logan(SheppLoganToft(); fovs=(200mm,200mm))


#=
From that collection we can compute images, sinograms and spectra.
Use `phantom` to make digital images.
It is convenient (but not required)
to use `ImageGeoms` to help with the sampling.
=#

ig = ImageGeom(dims=(200,256), deltas=(1mm,1mm))
image = phantom(axes(ig)..., objects)
jim(axes(ig), image; xlabel="x", ylabel="y", title="SheppLoganToft")


# Here is the sinogram corresponding to this phantom,
# computed analytically from the ellipse parameters using `radon`:

r = LinRange(-100mm,100mm,401)
ϕ = deg2rad.(0:180)
sino = radon(r, ϕ, objects)
jim(r, ϕ, sino; title="Sinogram", xlabel="r", ylabel="ϕ")


#=
## Spectra
Here is the 2D spectrum (Fourier transform) of this phantom,
computed analytically from the ellipse parameters using `spectrum`:
=#

kspace = spectrum(axesf(ig)..., objects)
jim(axesf(ig), log10.(abs.(kspace/(1mm)^2)); xlabel="ν₁", ylabel="ν₂", title="log10|Spectrum|")

#=
The 2D Fourier transform formula used here is:
```math
F(ν_1, ν_2) = ∫ ∫ f(x,y) \, \mathrm{e}^{-ı 2π (ν_1 x + ν_2 y)} \, \mathrm{d} x \, \mathrm{d} y,
```
so if ``x`` and ``y`` have units mm (for example),
then the units of the spatial frequency variables
``ν_1`` and ``ν_2`` are cycles/mm.
=#


#=
## 2D Rotation

All of the 2D objects (ellipses etc.)
in this package can be rotated by an angle ``ϕ``.

This package treats the rotation angle ``ϕ``
as defining a
[rotation of the object](https://en.wikipedia.org/wiki/Rotation_(mathematics)#Two_dimensions).
Be aware that a
[rotation of the axes](https://en.wikipedia.org/wiki/Rotation_of_axes)
has the opposite sign convention.

After rotation by ``ϕ``,
any point ``(x,y)`` in the original ellipse
becomes the point
```math
\left[ \begin{matrix}
x' \\ y'
\end{matrix} \right]
=
\left[ \begin{matrix}
\cos(ϕ) & -\sin(ϕ) \\ \sin(ϕ) & \cos(ϕ)
\end{matrix} \right]
\left[ \begin{matrix}
x \\ y
\end{matrix} \right],
```
as illustrated by the blue star below
when rotating an ellipse
by ``ϕ = π/6``.
=#

ellipse0 = ellipse(0, 0, 8, 4, 0, 1)
ϕ1s = :(π/6)
ϕ1 = eval(ϕ1s)
ellipse1 = ellipse(0, 0, 8, 4, ϕ1, 1)

x = LinRange(-9, 9, 181)
y = LinRange(-8, 8, 161)
pic0 = phantom(x, y, [ellipse0])
pic1 = phantom(x, y, [ellipse1])

p0 = jim(x, y, pic0, "Original ellipse"; prompt=:false)
x0,y0 = 7,0
scatter!([x0], [y0], marker=:star, color=:blue)
point1 = [cos(ϕ1) -sin(ϕ1); sin(ϕ1) cos(ϕ1)] * [x0; y0] # rotate point
x1,y1 = point1[1], point1[2]
p1 = jim(x, y, pic1, "Rotated by ϕ = $ϕ1s"; prompt=:false)
scatter!([x1], [y1], marker=:star, color=:blue)
jim(p0, p1)



# ## Reproducibility

# This page was generated with the following version of Julia:

io = IOBuffer(); versioninfo(io); split(String(take!(io)), '\n')


# And with the following package versions

import Pkg; Pkg.status()
