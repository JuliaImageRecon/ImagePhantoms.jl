#=
# [ImagePhantoms overview](@id 01-overview)

This page explains the Julia package
[`ImagePhantoms`](https://github.com/JuliaImageRecon/ImagePhantoms.jl).
=#

#srcURL


# ### Setup

# Packages needed here.

#nb import Pkg
#nb Pkg.add(["ImageGeoms", "MIRTjim", "Unitful"])

using ImagePhantoms
using ImageGeoms: ImageGeom, axesf
using MIRTjim: jim, prompt, mid3
using Plots; default(markerstrokecolor=:auto, label="")
using Unitful: mm
using InteractiveUtils: versioninfo

# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);


#=
## Overview

When developing image reconstruction methods,
it can be helpful to simulate data (e.g., sinograms)
using software-defined images called phantoms.

The simplest method here is to make a Shepp-Logan phantom image
similar its use in other packages.
=#

image = shepp_logan(256) # CT version by default
jim(image, "SheppLogan"; clim=(0.9, 1.1))


#=
## Sinograms

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

r = range(-100mm,100mm,401)
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

x = range(-9, 9, 181)
y = range(-8, 8, 161)
pic0 = phantom(x, y, [ellipse0])
pic1 = phantom(x, y, [ellipse1])

marker = :star
p0 = jim(x, y, pic0, "Original ellipse";
    xlabel="x", ylabel="y", size=(700,300), prompt=:false)
x0,y0 = 7,0
scatter!([x0], [y0], color=:blue; marker)
point1 = [cos(ϕ1) -sin(ϕ1); sin(ϕ1) cos(ϕ1)] * [x0; y0] # rotate point
x1,y1 = point1[1], point1[2]
p1 = jim(x, y, pic1, "Rotated by ϕ = $ϕ1s";
    xlabel="x", ylabel="y", size=(700,300), prompt=:false)
scatter!([x1], [y1], color=:blue; marker)
jim(p0, p1)


#=
## 3D Rotation

For a 3D object,
there are three rotation angles,
often called
[Euler angles](https://en.wikipedia.org/wiki/Euler_angles),
and there are many possible conventions
for the names and ordering of these angles.

This package denotes
the three angles as ``ϕ,θ,ψ.``

### Rotation about ``z`` by ``ϕ``

For consistency with the 2D case,
the first of the three angles,
``ϕ``,
denotes rotation in the ``(x,y)`` plane,
i.e., around the ``z``-axis.
In
[wikipedia's notation](https://en.wikipedia.org/w/index.php?title=Rotation_matrix&section=9#In_three_dimensions)
this is
``R_z(ϕ)``,
defined as
```math
\left[ \begin{matrix}
x' \\ y' \\ z'
\end{matrix} \right]
=
\left[ \begin{matrix}
\cos(θ) & -\sin(θ) & 0
\\
\sin(θ) & \cos(θ) & 0
\\
0 & 0 & 1
\\
\end{matrix} \right]
\left[ \begin{matrix}
x \\ y \\ z
\end{matrix} \right]
```
as illustrated by the blue star below.

Here is an illustration
for $ϕ = π/6$.
=#

ellipsoid0 = ellipsoid((0, 0, 0), (8, 4, 2), (0, 0, 0), 1)
ϕ1s = :(π/6)
ϕ1 = eval(ϕ1s)
ellipsoid1 = ellipsoid((0, 0, 0), (8, 4, 2), (ϕ1, 0, 0), 1)

x = range(-9, 9, 181)
y = range(-8, 8, 161)
z = [0]
pic0 = phantom(x, y, z, [ellipsoid0])
pic1 = phantom(x, y, z, [ellipsoid1])

p0z = jim(x, y, pic0,
    "Original ellipsoid:\n(x,y) slice";
    xlabel="x", ylabel="y", size=(700,300), prompt=:false)
x0,y0 = 7,0
scatter!([x0], [y0], color=:blue; marker)
Rz(ϕ) = [cos(ϕ) -sin(ϕ) 0; sin(ϕ) cos(ϕ) 0; 0 0 1]
point1 = Rz(ϕ1) * [x0; y0; 0] # rotate point
x1,y1 = point1[1], point1[2]
p1z = jim(x, y, pic1,
    "Rotated about z\nby ϕ = $ϕ1s\n(z out of 'board')";
    xlabel="x", ylabel="y", size=(700,350), prompt=:false)
scatter!([x1], [y1], color=:blue; marker)
jim(p0z, p1z)


#=
### Rotation about ``y`` by ``θ``

The 2nd of the three angles,
``θ``,
corresponds to rotation around the ``y``-axis,
which has the opposite sign
when using the
[right hand rule](https://en.wikipedia.org/wiki/Right-hand_rule#A_rotating_body):
```math
\left[ \begin{matrix}
x' \\ y' \\ z'
\end{matrix} \right]
=
\left[ \begin{matrix}
\cos(θ) & 0 & \sin(θ)
\\
0 & 1 & 0
\\
-\sin(θ) & 0 & \cos(θ)
\\
\end{matrix} \right]
\left[ \begin{matrix}
x \\ y \\ z
\end{matrix} \right],
```
as illustrated by the green star below.

In
[wikipedia notation](https://en.wikipedia.org/w/index.php?title=Rotation_matrix&section=9#In_three_dimensions)
this is
``R_y(θ)``.

Here is an illustration
of rotating an ellipsoid
for ``θ = π/6``.
=#

ellipsoid0 = ellipsoid((0, 0, 0), (8, 4, 2), (0, 0, 0), 1)
θ1s = :(π/6)
θ1 = eval(θ1s)
ellipsoid1 = ellipsoid((0, 0, 0), (8, 4, 2), (0, θ1, 0), 1)

x = range(-9, 9, 181)
y = [0]
z = range(-8, 8, 161)
pic0 = phantom(x, y, z, [ellipsoid0]); pic0 = selectdim(pic0, 2, 1)
pic1 = phantom(x, y, z, [ellipsoid1]); pic1 = selectdim(pic1, 2, 1)

p0y = jim(x, z, pic0,
    "Original ellipsoid:\n (x,z) slice";
    xlabel="x", ylabel="z", size=(700,350), prompt=false)
x0,z0 = 7,0
scatter!([x0], [z0], color=:green; marker)
Ry(θ) = [cos(θ) 0 sin(θ); 0 1 0; -sin(θ) 0 cos(θ)]
point1 = Ry(θ1) * [x0; 0; z0] # rotate point
x1,z1 = point1[1], point1[3]
p1y = jim(x, z, pic1,
    "Rotated about y\nby θ = $θ1s\n(y into 'board')";
    xlabel="x", ylabel="z", size=(700,350), prompt=false)
scatter!([x1], [z1], color=:green; marker)
jim(p0y, p1y)


#=
### Rotation about ``x`` by ``ψ``

The 3rd of the three angles,
``ψ``,
corresponds to rotation around the ``x``-axis:

```math
\left[ \begin{matrix}
x' \\ y' \\ z'
\end{matrix} \right]
=
\left[ \begin{matrix}
1 & 0 & 0
\\
0 & \cos(ψ) & -\sin(ψ)
\\
0 & \sin(ψ) & \cos(ψ)
\\
\end{matrix} \right]
\left[ \begin{matrix}
x \\ y \\ z
\end{matrix} \right],
```
as illustrated by the red star below.

In
[wikipedia notation](https://en.wikipedia.org/w/index.php?title=Rotation_matrix&section=9#In_three_dimensions)
this is
``R_x(ψ)``.

Here is an illustration
of rotating an ellipsoid
for ``ψ = π/6``.
=#

ellipsoid0 = ellipsoid((0, 0, 0), (8, 4, 2), (0, 0, 0), 1)
ψ1s = :(π/6)
ψ1 = eval(ψ1s)
ellipsoid1 = ellipsoid((0, 0, 0), (8, 4, 2), (0, 0, ψ1), 1)

x = [0]
y = range(-9, 9, 181)
z = range(-8, 8, 161)
pic0 = phantom(x, y, z, [ellipsoid0]); pic0 = selectdim(pic0, 1, 1)
pic1 = phantom(x, y, z, [ellipsoid1]); pic1 = selectdim(pic1, 1, 1)

p0x = jim(y, z, pic0,
    "Original ellipsoid:\n (y,z) slice)";
    xlabel="y", ylabel="z", size=(700,350), prompt=false)
y0,z0 = 3,0
scatter!([y0], [z0], color=:red; marker)
Rx(ψ) = [1 0 0 ; 0 cos(ψ) -sin(ψ); 0 sin(ψ) cos(ψ)]
point1 = Rx(ψ1) * [0; y0; z0] # rotate point
y1,z1 = point1[2], point1[3]
p1x = jim(y, z, pic1,
    "Rotated about x\nby ψ = $ψ1s\n(x out of 'board')";
    xlabel="y", ylabel="z", size=(700,350), prompt=false)
scatter!([y1], [z1], color=:red; marker)
jim(p0x, p1x)


#src # Summarizing all in one plot:
#src jim(p0x, p0y, p0z, p1x, p1y, p1z)


#=
The remaining issue is the multiplication order
for multiple rotations.

To address that,
we first describe
how phantoms are generated in this package.
Every phantom shape
starts with a base function
that is translated, rotated, and scaled
to make the final object.
For example,
an ellipsoid is a transformed sphere.
Specifically,
if ``e(r)`` is the ellipsoid function,
and ``s(r)`` is the unit sphere function,
where
``r = (x,y,z)``,
then
```math
e(r) = s( (R_{xyz}^T (r - c)) ⊘ w )
```
where ``c`` is the `center` of the ellipsoid,
and ``r`` is the `width` parameter (actually radii),
``⊘`` denotes element-wise division,
and ``R_{xyz}^T``
is the *inverse*
of the 3D rotation matrix
``R_{xyz}``
defined by
``R_{xyz} = R_x(ψ) R_y(θ) R_z(ϕ)``.

Note that ``r`` and ``c`` and ``w``
all must have identical units,
so the argument
``(R_{xyz}^T (r - c)) ⊘ w``
passed the unit-sphere function
is unitless,
as it must be.
The non-exported `phantom1` function
for each shape
defines the base shape function.
For the unit sphere it is simply
`sum(abs2, r) ≤ 1`.

Rearranging the equation
``r' = (R_{xyz}^T (r - c)) ⊘ w``
yields
``r = R_{xyz} (w ⊙ r') + c``.
So the process of transforming a unit sphere
to an ellipsoid starts with scaling,
then rotation by
``R_{xyz}``,
which rotates first around the ``z``-axis,
and then finally translating.

Here is an illustration
where one can see that all three axes were rotated.
=#

ellipsoid0 = ellipsoid((0, 0, 0), (8, 4, 2), (0, 0, 0), 1)
ϕ1s = :(π/6)
ϕ1 = eval(ϕ1s)
θ1s = :(π/7)
θ1 = eval(θ1s)
ψ1s = :(π/8)
ψ1 = eval(ψ1s)
ellipsoid1 = ellipsoid((0, 0, 0), (8, 4, 2), (ϕ1, θ1, ψ1), 1)

x = range(-9, 9, 181)
y = range(-8, 8, 161)
z = range(-7, 7, 71)
pic0 = phantom(x, y, z, [ellipsoid0])
pic1 = phantom(x, y, z, [ellipsoid1])

p0a = jim(mid3(pic0),
    "Original ellipsoid\n(central slices)";
    xlabel="x", ylabel="y", size=(700,320), prompt=:false)
p1a = jim(mid3(pic1),
    "Rotated\nϕ = $ϕ1s, θ = $θ1s, ψ = $ψ1s";
    xlabel="x", ylabel="y", size=(700,320), prompt=:false)
jim(p0a, p1a)


#=
## Spectra rotation

The `spectrum` method
accounts for the translation, rotation, and scaling
of the base shape function
using elementary Fourier transform properties.

The following code
first shows the spectra of an ellipsoid
before and after rotating it.
=#

ig = ImageGeom(length.((x,y,z)), map(x -> x[2]-x[1], (x,y,z)))
kspace0 = spectrum(axesf(ig)..., [ellipsoid0])
kspace1 = spectrum(axesf(ig)..., [ellipsoid1])
@assert ImagePhantoms.volume(ellipsoid0) == ImagePhantoms.volume(ellipsoid1)

clim = (-6, 0)
p0s = jim(axesf(ig), log10.(abs.(kspace0 / ImagePhantoms.volume(ellipsoid0)));
    clim, xlabel="ν₁", ylabel="ν₂", title="log10|Spectrum original ellipsoid|")

#
p1s = jim(axesf(ig), log10.(abs.(kspace1 / ImagePhantoms.volume(ellipsoid1)));
    clim, xlabel="ν₁", ylabel="ν₂", title="log10|Spectrum rotated ellipsoid|")

# The following code verifies that the rotated spectrum matches
k = Iterators.product(axesf(ig)...) # tuples of (kx,ky,kz) values
R = Rx(ψ1) * Ry(θ1) * Rz(ϕ1) # Rxyz rotation matrix
kr = (tuple((R' * collect(k))...) for k in k) # rotate each k-space tuple
kspace0r = spectrum(kr, [ellipsoid0]) # evaluate spectrum at rotated tuples
@assert kspace0r ≈ kspace1


include("../../../inc/reproduce.jl")
