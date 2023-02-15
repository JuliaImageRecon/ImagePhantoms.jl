#=
# [3D Shepp-Logan Phantom](@id 37-shepp3)

This page illustrates the 3D Shepp-Logan phantom(s)
in the Julia package
[`ImagePhantoms`](https://github.com/JuliaImageRecon/ImagePhantoms.jl).
=#

#srcURL


# ### Setup

# Packages needed here.

using ImagePhantoms: ellipsoid_parameters, ellipsoid, phantom
using ImageGeoms: ImageGeom, axes
using MIRTjim: jim, prompt
using Unitful: g, cm


# The following line is helpful when running this file as a script;
# it prompts user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);


#=
## Overview

Currently this package provides
one version of a 3D Shepp-Logan phantom
consisting of ellipsoids.

For completeness,
we illustrate it using units,
though units are not required.

We use cm for the spatial units
and g/cm³ (density)
for the tissue values.
=#

fovs = (24cm, 24cm, 20cm)
u = (1, 1, 1*g/cm^3)
params = ellipsoid_parameters( ; fovs, u)
ob = ellipsoid(params); # Vector of Ellipsoid objects


#=
To visualize this phantom,
we sample it
with the help of `ImageGeoms`,
using over-sampling
to account for partial volume effects.
=#

dims = (128,130,30)
ig = ImageGeom( ; dims, deltas = fovs ./ dims )

oversample = 3
image = phantom(axes(ig)..., ob, oversample)
clim = (0.95, 1.05)
jim(axes(ig)[1:2]..., image; title = "3D Shepp-Logan phantom slices", clim)
