# ImagePhantoms.jl

https://github.com/JuliaImageRecon/ImagePhantoms.jl (_WIP!_)

[![action status][action-img]][action-url]
[![pkgeval status][pkgeval-img]][pkgeval-url]
[![codecov][codecov-img]][codecov-url]
[![license][license-img]][license-url]
[![docs-stable][docs-stable-img]][docs-stable-url]
[![docs-dev][docs-dev-img]][docs-dev-url]
[![code-style][code-blue-img]][code-blue-url]

This Julia language repo
provides tools for working with software-defined image phantoms
like the Shepp-Logan phantom.

For explanations see the documentation
using the blue "docs" links above.


### Getting started

```julia
using Pkg
Pkg.add("ImagePhantoms")
```


### Example

```julia
using ImagePhantoms
using ImageGeoms # not required, but helpful
todo
```

### Documentation

For more examples with graphics,
see the
[documentation][docs-stable-url].


### Dependents

* [Michigan Image Reconstruction Toolbox (MIRT)](https://github.com/JeffFessler/MIRT.jl)
* [Sinograms.jl](https://github.com/todo/Sinograms.jl)
* [SPECTrecon.jl](https://github.com/todo/SPECTrecon.jl)
* See [juliahub](https://juliahub.com/ui/Search?q=ImagePhantoms&type=packages)


### Related packages

* [AxisArrays](https://github.com/JuliaArrays/AxisArrays.jl)
* [ImageGeoms](https://github.com/JuliaImageRecon/ImageGeoms.jl)
* [JuliaImages/Images.jl](https://github.com/JuliaImages/Images.jl) `shepp_logan`
* [TestImages.jl](https://github.com/JuliaImages/TestImages.jl): `shepp_logan`


### Compatibility

Tested with Julia ≥ 1.6.

<!-- URLs -->
[action-img]: https://github.com/JuliaImageRecon/ImagePhantoms.jl/workflows/CI/badge.svg
[action-url]: https://github.com/JuliaImageRecon/ImagePhantoms.jl/actions
[build-img]: https://github.com/JuliaImageRecon/ImagePhantoms.jl/workflows/CI/badge.svg?branch=main
[build-url]: https://github.com/JuliaImageRecon/ImagePhantoms.jl/actions?query=workflow%3ACI+branch%3Amain
[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/I/ImagePhantoms.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/I/ImagePhantoms.html
[code-blue-img]: https://img.shields.io/badge/code%20style-blue-4495d1.svg
[code-blue-url]: https://github.com/invenia/BlueStyle
[codecov-img]: https://codecov.io/github/JuliaImageRecon/ImagePhantoms.jl/coverage.svg?branch=main
[codecov-url]: https://codecov.io/github/JuliaImageRecon/ImagePhantoms.jl?branch=main
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://JuliaImageRecon.github.io/ImagePhantoms.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://JuliaImageRecon.github.io/ImagePhantoms.jl/dev
[license-img]: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat
[license-url]: LICENSE
