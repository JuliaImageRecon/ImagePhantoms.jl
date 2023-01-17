```@meta
CurrentModule = ImagePhantoms
```

# ImagePhantoms.jl Documentation

## Overview

This Julia package
[`ImagePhantoms`](https://github.com/JuliaImageRecon/ImagePhantoms.jl)
provides tools
for constructing digital software "phantoms"
used for testing image reconstruction algorithms.
The most famous such phantom
is the "Shepp Logan" phantom
from this
[1974 paper](https://doi.org/10.1109/TNS.1974.6499235).
(A variant of that phantom serves
as the
logo for the
[JuliaImageRecon](https://github.com/JuliaImageRecon)
suite of tools.)

A phantom is a collection (iterable) of shapes
(e.g., ellipses, rectangles).
One purpose of this package
is to avoid
the
[inverse crime](https://doi.org/10.1016/j.cam.2005.09.027)
of using a discretized or sampled image
to generate sinograms or spectra.

The shapes in this package
have methods useful
for simulating data:
* `phantom` returns a function of `(x,y)` or `(x,y,z)`
   that one can evaluate on a grid to make a pixelated or voxelized phantom.
* `radon` returns a function of `(r,ϕ)`
   that one can evaluate on a grid
   to make a sampled 2D parallel-beam sinogram,
   or evaluate appropriately to make a sampled fan-beam sinogram.
   For 3D objects,
  `radon` returns a function of `(u,v,ϕ,θ)`
   that one can evaluate to compute projection views.
   See the "3D geometry" example
   for description of the coordinate system.
* `spectrum` returns a function of spatial frequencies
  that one can sample to simulate k-space data (e.g., in MRI).

See the Examples tab for details.

The
[Michigan Image Reconstruction Toolbox (MIRT)](https://github.com/JeffFessler/MIRT.jl)
currently has an older interface `ellipse_im`, `rect_im`, etc.,
similar to the functions of the same name in the
[Matlab version of MIRT](https://github.com/JeffFessler/mirt)
provided for backward compatibility.
Using `ImagePhantoms` is recommended for Julia work.


## Units

This package allows the shapes to be described
with physical units,
e.g., using
[Unitful.jl](https://github.com/PainterQubits/Unitful.jl).
Using units is recommended but not required.
If units are not used explicitly,
then the user
must be especially careful
to use values
corresponding to consistent units.
For example,
if the shape sizes are in cm,
then
* the arguments to `phantom`
  must also be in cm units,
* the first (`r`) argument to `radon`
  must also be in cm
  and the `ϕ` argument must be in radians,
* the spatial frequency arguments to `spectrum`
  must be in cycles/cm units.
