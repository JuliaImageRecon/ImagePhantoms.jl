#=
# [3D geometry](@id 30-3d)

This page explains the 3D X-ray transform geometry
for the Julia package
[`ImagePhantoms`](https://github.com/JuliaImageRecon/ImagePhantoms.jl).
=#

#srcURL


# ## Setup

# Packages needed here.

#src using ImagePhantoms
#src using ImageGeoms: ImageGeom, axesf
#src using MIRTjim: jim, prompt
#src using Unitful: mm


#src # The following line is helpful when running this file as a script;
#src # this way it will prompt user to hit a key after each figure is displayed.

#src isinteractive() ? jim(:prompt, true) : prompt(:draw);


#=
## Overview

Most practical tomographic imaging problems are 3D,
so this package also supports calculating
3D line integrals through 3D phantom objects.


## 3D projection geometry

Given a 3D phantom object,
the function `radon`
returns a 4-argument function
having arguments
``(u,v,ϕ,θ)``,
where ``u,v`` denote the coordinates
on a 2D projection plane,
``ϕ`` denotes the
[azimuthal angle](https://en.wikipedia.org/wiki/Azimuth)
angle
and
``θ`` denotes the
[polar angle](https://en.wikipedia.org/wiki/Spherical_coordinate_system).

The mathematical definition of the
[Radon transform of a 3D function](https://en.wikipedia.org/wiki/Radon_transform)
is a collection of (2D) plane integrals,
whereas the
[X-ray transform](https://en.wikipedia.org/wiki/X-ray_transform)
is a collection of (1D) line integrals.
See
[Section II.1, Natterer 2001](https://doi.org/10.1137/1.9780898719284).
So strictly speaking the `radon` function is a misnomer in 3D,
whereas for 2D functions
the Radon transform and the X-ray transform coincide.

The coordinate system used here is defined as follows.
Start by defining a point on the "detector plane" as follows:
```math
\vec{p} = \vec{p}(u, v, ϕ, θ)
= u \vec{e}_1 + v \vec{e}_3
= (u \cos ϕ + v \sin ϕ \sin θ, u \sin ϕ - v \cos ϕ \sin θ, v \cos θ) ∈ ℝ^3,
```
where
```math
\vec{e}_1 = (\cos ϕ, \sin ϕ, 0)
,\qquad
\vec{e}_3 = (\sin ϕ \sin θ, -\cos ϕ \sin θ, \cos θ).
```

Now define the (X-ray) projection
of a 3D object ``f(\vec{x})``
as
```math
p(u, v, ϕ, θ)
= ∫ f(\vec{p} + ℓ \, \vec{e}) \, \mathrm{d} ℓ
= ∫_{-∞}^{∞} f( \vec{p}(u,v,ϕ,θ) + ℓ \, \vec{e}(ϕ,θ) ) \, \mathrm{d} ℓ,
```
where
```math
\vec{e}(ϕ,θ) = (-\sin ϕ \cos θ, \cos ϕ \cos θ, \sin θ).
```

When ``θ=0``,
then the X-ray transform
``p(u, v, ϕ, θ)``
is a collection of 2D sinograms,
one for each slice of ``f(x,y,z)``:
```math
p(u, v, ϕ, 0)
= ∫_{-∞}^{∞} f(u \cos ϕ - ℓ \sin ϕ, u \sin ϕ + ℓ \cos ϕ, v) \, \mathrm{d} ℓ.
```

For example projection views,
see the Ellipsoid examples.

=#
