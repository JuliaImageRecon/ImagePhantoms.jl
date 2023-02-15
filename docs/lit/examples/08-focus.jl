#=
# [Focus Chart](@id 08-focus)

This page illustrates the `focus_chart` method in the Julia package
[`ImagePhantoms`](https://github.com/JuliaImageRecon/ImagePhantoms.jl).
=#

#srcURL


# ### Setup

# Packages needed here.

using ImagePhantoms: focus_chart, phantom
using MIRTjim: jim


#=
## Focus chart phantom image

One of the shapes in `ImagePhantoms` is an equilateral triangle,
and by scaling and translating this shape
one can define a focus chart phantom.
=#

ob = focus_chart( ; nspoke = 56, value = 4)
x = range(-1,1,2^9) * 1.1
y = x
image = phantom(ob).(x, y')
jim(x, y, image; title = "Focus chart phantom")
