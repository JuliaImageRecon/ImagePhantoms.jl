```@meta
CurrentModule = ImagePhantoms
```

# ImagePhantoms.jl Documentation

## Overview

This Julia package provides tools
for constructing digital software "phantoms"
used for testing image reconstruction algorithms.
The most famous such phantom
is the "Shepp Logan" phantom
from this
[1974 paper](http://doi.org/10.1109/TNS.1974.6499235).
(A variant of that phantom serves
as the
logo for the
[JuliaImageRecon](https://github.com/JuliaImageRecon)
suite of tools.

A phantom is a collection (iterable) of shapes
(e.g., ellipses, rectangles).
This package allows the shapes to be described
with physical units,
e.g., using
[Unitful.jl](https://github.com/PainterQubits/Unitful.jl).

The shapes in this package
have methods useful
for simulating data:
* `phantom` returns a function of `(x,y)` or `(x,y,z)`
   that one can evaluate on a grid to make a pixelated or voxelized phantom.
* `radon` returns a function of `(r,ϕ)`
   that one can evaluate on a grid
   to make a sampled 2D parallel-beam sinogram,
   or evaluate appropriately to make a sampled fan-beam sinogram.
* `spectrum` returns a function of spatial frequencies
  that one can sample to simulate k-space data (e.g., in MRI).

See the
[Examples](@ref 1-overview)
tab to the left for details.

The
[Michigan Image Reconstruction Toolbox (MIRT)](https://github.com/JeffFessler/MIRT.jl)
currently has an older interface `ellipse_im`, `rect_im`, etc.,
similar to the functions of the same name in the
[Matlab version of MIRT](https://github.com/JeffFessler/mirt)
provided for backward compatibility.
Using `ImagePhantoms` is recommended for Julia work.
