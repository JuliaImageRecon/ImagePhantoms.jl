#---------------------------------------------------------
# # [ImagePhantoms overview](@id 1-overview)
#---------------------------------------------------------

# This page explains the Julia package
# [`ImagePhantoms`](https://github.com/JuliaImageRecon/ImagePhantoms.jl).

# ### Setup

# Packages needed here.

using ImagePhantoms
using ImageGeoms: ImageGeom
using MIRTjim: jim, prompt
using Plots: scatter, plot!, default; default(markerstrokecolor=:auto)

# The following line is helpful when running this example.jl file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);

# ### Overview

# When developing image reconstruction methods,
# it can be helpful to simulate data (e.g., sinograms)
# using software-defined images called phantoms.

Ellipse(todo)

# The key parameters of a grid of image pixels are
# * the size (dimensions) of the grid, e.g., `128 × 128`,
# * the spacing of the pixels, e.g., `1mm × 1mm`,
# * the offset of the pixels relative to the origin, e.g., `(0,0)`

# The data type `ImageGeom` describes such a geometry,
# for arbitrary dimensions (2D, 3D, etc.).

# There are several ways to construct this structure.
# The default is a `128 × 128` grid
# with pixel size ``\Delta_X = \Delta_Y = 1`` (unitless) and zero offset:

ig = ImageGeom()


# Here is a 3D example with non-cubic voxel size:

ig = ImageGeom( (512,512,128), (1,1,2), (0,0,0) )


# To avoid remembering the order of the arguments,
# named keyword pairs are also supported:

ig = ImageGeom( dims=(512,512,128), deltas=(1,1,2), offsets=(0,0,0) )


# ### Units

# The pixel dimensions `deltas` can (and should!) be values with units.

# Here is an example for a video (2D+time) with 12 frames per second:
using UnitfulRecipes
using Unitful: mm, s

ig = ImageGeom( dims=(640,480,1000), deltas=(1mm,1mm,(1//12)s) )


# ### Methods

# An ImageGeom object has quite a few methods;
# `axes` and `axis` are especially useful:

ig = ImageGeom( dims=(7,8), deltas=(3,2), offsets=(0,0.5) )
axis(ig, 2)


# For an axis of length `n` with spacing `Δ` (possibly with units)
# and (always unitless but possibly non-integer) `offset` the axis
# is a subtype of `AbstractRange` of the form
# `( (0:n-1) .- ((n - 1)/2 + offset) ) * Δ`


# These axes are useful for plotting:
ig = ImageGeom( dims=(12,10), deltas=(1mm,1mm), offsets=(0.5,0.5) )


#
showgrid = (ig) -> begin # x,y grid locations of pixel centers
    x = axis(ig, 1)
    y = axis(ig, 2)
	(xg, yg) = grids(ig)
    scatter(xg, yg, label="", xlabel="x", ylabel="y",
        xlims = maximum(abs, x) * 1.2 .* (-1,1),
        xticks = [x[1], zero(eltype(x)), x[end]],
        ylims = maximum(abs, y) * 1.2 .* (-1,1),
        yticks = [y[1], zero(eltype(y)), y[end]],
        aspect_ratio=1, title="offsets $(ig.offsets)")
end
showgrid(ig)


# Unit labels on the axes due to `UnitfulRecipes.jl`
prompt();


# ### Offsets (unitless translation of grid)

# The default `offsets` are zeros,
# corresponding to symmetric sampling around origin:

ig = ImageGeom( dims=(12,10), deltas=(1mm,1mm) )
p = showgrid(ig)

#
prompt();


# That default for `offsets` is natural for tomography
# when considering finite pixel size:

square = (x,y,Δ) -> plot!(p, label="", color=:black,
    x .+ Δ[1] * ([0,1,1,0,0] .- 0.5),
    y .+ Δ[2] * ([0,0,1,1,0] .- 0.5),
)
square2 = (x,y) -> square(x, y, ig.deltas)
square2.(grids(ig)...)
plot!(p)


#
prompt();


# In that default geometry, the center `(0,0)` of the image
# is at a corner of the middle 4 pixels (for even image sizes).
# That default is typical for tomographic imaging (e.g., CT, PET, SPECT).
# One must be careful when using operations like `imrotate` or `fft`.



# ### AxisArrays

# There is a natural connection between `ImageGeom` and `AxisArrays`.
# Note the automatic labeling of units (when relevant) on all axes by
# [MIRTjim.jim](https://github.com/JeffFessler/MIRTjim.jl).

using AxisArrays
using Unitful: mm
ig = ImageGeom( dims=(60,48), deltas=(1.5mm,1mm) )
za = AxisArray( ellipse(ig) * 10/mm ; x=axis(ig,1), y=axis(ig,2) )
jim(za, "AxisArray example")

#
prompt();


# ### Resizing

# Often we have a target grid in mind but want coarser sampling for debugging.
# The `downsample` method is useful for this.

ig = ImageGeom( dims = (512,512), deltas = (500mm,500mm) ./ 512 )
ig_down = downsample(ig, 4)


# Other times we want to avoid an "inverse crime" by using finer sampling
# to simulate data; use `oversample` for this.

ig_over = oversample(ig, (2,2))
