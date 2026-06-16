# ImagePhantoms.jl

https://github.com/JuliaImageRecon/ImagePhantoms.jl
<img src="https://github.com/JuliaImageRecon/ImagePhantoms.jl/blob/main/docs/src/assets/logo.svg" alt="logo" width="150"/>

[![docs-stable][docs-stable-img]][docs-stable-url]
[![docs-dev][docs-dev-img]][docs-dev-url]
[![action][action-img]][action-url]
[![Aqua QA][aqua-img]][aqua-url]
[![codecov][codecov-img]][codecov-url]
[![deps][deps-img]][deps-url]
[![license][license-img]][license-url]
[![pkgeval][pkgeval-img]][pkgeval-url]
[![version][ver-img]][ver-url]

This Julia language repo
provides tools for working with software-defined image phantoms
like the Shepp-Logan phantom
(both 2D and 3D versions).


For each phantom shape,
methods are available for computing samples of
its image
(using `phantom`),
its Radon transform (2D) or X-ray line integrals (3D),
(using `radon`),
and its 2D or 3D Fourier transform
(using `spectrum`).

For explanations and examples
see the documentation
using the blue "docs" links above.


### Getting started

```julia
using Pkg
Pkg.add("ImagePhantoms")
```


### Example

```julia
using ImagePhantoms
image = shepp_logan(256)

p = shepp_logan(SheppLoganToft())
sinogram = radon(range(-1,1,101), deg2rad.(0:180), p) # 101 × 181
```


### Example

```julia
using ImagePhantoms
using MIRTjim: jim
using Plots # @animate, gif
function disk_phantom(title::String)
    (dx,dy) = (1, 1)
    (M,N) = (2^8,2^8)
    x = (-M÷2:M÷2-1) * dx
    y = (-N÷2:N÷2-1) * dy
    params = disk_phantom_params( ; rhead = () -> rand(100:105))
    objects = ellipse(params) # vector of Object2d{Ellipse}
    img = phantom(x, y, objects) # sampled at all (x,y) pairs
    jim(x, y, img; title, clim=(0,1300))
end
anim = @animate for i in 1:8
    disk_phantom("Realization $i")
end
gif(anim, "disk.gif", fps = 8)
```

![animated phantom gif](https://github.com/JuliaImageRecon/ImagePhantoms.jl/blob/gh-pages/dev/generated/examples/disk.gif)


### Philosophy

Often "phantoms" are treated as digital images.
Here, the shapes (rectangles, gaussians, ellipses, etc.)
are all defined *analytically*
as functions,
as are their line integrals
and Fourier transforms.
Then one can sample those analytical functions
to make digital images, sinograms, and spectra.


### Parallel MRI (SENSE)

Most of the methods here are of general use
for any imaging modality.
There are a few methods
(`mri_smap_basis`, `mri_smap_fit`, `mri_spectra`)
that are specific to simulating parallel MRI
(multiple receive coils)
based on the
[2012 paper by Guerquin-Kern et al.](https://doi.org/10.1109/TMI.2011.2174158).
See the
[documentation][docs-stable-url]
for details.


### Documentation

For more examples with graphics,
see the
[documentation][docs-stable-url].


Currently the package supports
the following 2D shapes:
ellipses/circles, rectangles/squares, gaussians, triangles,
and the following 3D shapes:
ellipsoids/spheres, cuboids/cubes, gaussians, cylinders, cones.


### Dependents

* [Michigan Image Reconstruction Toolbox (MIRT)](https://github.com/JeffFessler/MIRT.jl)
* [Sinograms.jl](https://github.com/JuliaImageRecon/Sinograms.jl)
* [SPECTrecon.jl](https://github.com/JuliaImageRecon/SPECTrecon.jl)
* See [juliahub](https://juliahub.com/ui/Search?q=ImagePhantoms&type=packages)


### Related packages

* [AxisArrays.jl](https://github.com/JuliaArrays/AxisArrays.jl)
* [ImageGeoms.jl](https://github.com/JuliaImageRecon/ImageGeoms.jl)
* [JuliaImages/Images.jl](https://github.com/JuliaImages/Images.jl) `shepp_logan`
* [TestImages.jl](https://github.com/JuliaImages/TestImages.jl): `shepp_logan`
* [KomaMRI.jl](https://github.com/cncastillo/KomaMRI.jl)


### Compatibility

Tested with Julia ≥ 1.12.

<!-- URLs -->
[action-img]: https://github.com/JuliaImageRecon/ImagePhantoms.jl/workflows/CI/badge.svg
[action-url]: https://github.com/JuliaImageRecon/ImagePhantoms.jl/actions

[aqua-img]: https://juliatesting.github.io/Aqua.jl/dev/assets/badge.svg
[aqua-url]: https://github.com/JuliaTesting/Aqua.jl

[codecov-img]: https://codecov.io/github/JuliaImageRecon/ImagePhantoms.jl/coverage.svg
[codecov-url]: https://codecov.io/github/JuliaImageRecon/ImagePhantoms.jl

[deps-img]: https://juliahub.com/docs/ImagePhantoms/deps.svg
[deps-url]: https://juliahub.com/ui/Packages/ImagePhantoms

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://JuliaImageRecon.github.io/ImagePhantoms.jl/dev
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://JuliaImageRecon.github.io/ImagePhantoms.jl/stable

[license-img]: https://img.shields.io/badge/license-MIT-brightgreen.svg
[license-url]: LICENSE

[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/I/ImagePhantoms.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/I/ImagePhantoms.html

[ver-img]: https://juliahub.com/docs/ImagePhantoms/version.svg
[ver-url]: https://juliahub.com/ui/Packages/ImagePhantoms
