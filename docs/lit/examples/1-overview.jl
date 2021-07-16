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

# WIP - more later about Shepp-Logan phantom, with units

using UnitfulRecipes
using Unitful: mm
