#---------------------------------------------------------
# # [Focus Chart](@id 08-focus)
#---------------------------------------------------------

#=
This page illustrates the `focus_chart` method in the Julia package
[`ImagePhantoms`](https://github.com/JuliaImageRecon/ImagePhantoms.jl).

This page was generated from a single Julia file:
[08-focus.jl](@__REPO_ROOT_URL__/08-focus.jl).
=#

#md # In any such Julia documentation,
#md # you can access the source code
#md # using the "Edit on GitHub" link in the top right.

#md # The corresponding notebook can be viewed in
#md # [nbviewer](http://nbviewer.jupyter.org/) here:
#md # [`08-focus.ipynb`](@__NBVIEWER_ROOT_URL__/08-focus.ipynb),
#md # and opened in [binder](https://mybinder.org/) here:
#md # [`08-focus.ipynb`](@__BINDER_ROOT_URL__/08-focus.ipynb).


# ### Setup

# Packages needed here.

using ImagePhantoms: focus_chart, phantom
using MIRTjim: jim


# ### Focus chart phantom image

#=
One of the shapes in `ImagePhantoms` is an equilateral triangle,
and by scaling and translating this shape
one can define a focus chart phantom.
=#

ob = focus_chart( ; nspoke = 56, value = 4)
x = LinRange(-1,1,2^9) * 1.1
y = x
image = phantom(ob).(x, y')
jim(x, y, image; title = "Focus chart phantom")
