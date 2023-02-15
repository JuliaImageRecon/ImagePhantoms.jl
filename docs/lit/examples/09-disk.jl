#=
# [Random Disks](@id 09-disk)

This page illustrates the `disk_phantom_params` method in the Julia package
[`ImagePhantoms`](https://github.com/JuliaImageRecon/ImagePhantoms.jl).
=#

#srcURL


# ### Setup

# Packages needed here.

using ImagePhantoms: ellipse, phantom, disk_phantom_params
using ImageGeoms: ImageGeom
using MIRTjim: jim, prompt
using Plots # @animate, gif

# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);

#=
## Overview

For training machine-learning methods for image reconstruction,
it can be helpful to have a way to generate a family of phantoms
drawn from some common distribution,
especially for debugging
or when real ground-truth data is limited.
The `disk_phantom_params` function is one example of how
one can use the tools in this package to generate such phantoms.
=#


# ## A single disk phantom image

# Make a digital image of it using `phantom` and display it.
function disk_phantom(title::String)
    (dx,dy) = (1, 1)
    (M,N) = (2^8,2^8)
    x = (-M÷2:M÷2-1) * dx
    y = (-N÷2:N÷2-1) * dy
    params = disk_phantom_params( ; rhead = () -> rand(100:105))
    objects = ellipse(params) # vector of Object{Ellipse}
    oversample = 3
    img = phantom(x, y, objects, oversample)
    jim(x, y, img; title, clim=(0,1300))
end
disk_phantom("A single disk phantom realization")


# ## Several realizations

anim = @animate for i in 1:8
    disk_phantom("Realization $i")
end
gif(anim, "disk.gif", fps = 6)
