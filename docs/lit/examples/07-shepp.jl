#=
# [Shepp-Logan Phantoms](@id 07-shepp)

This page illustrates the Shepp-Logan phantoms in the Julia package
[`ImagePhantoms`](https://github.com/JuliaImageRecon/ImagePhantoms.jl).
=#

#srcURL


# ### Setup

# Packages needed here.

using ImagePhantoms: shepp_logan, SouthPark
using ImagePhantoms: SheppLoganToft, SheppLoganEmis, SheppLoganBrainWeb
using ImagePhantoms: ellipse_parameters, ellipse, phantom
using MIRTjim: jim, prompt
using Unitful: g, cm


# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);


# ## Overview

# There are several variations of the Shepp-Logan phantom available.


#=
## CT version

Original version from:
Larry A Shepp, Benjamin F Logan,
"The Fourier reconstruction of a head section,"
IEEE Transactions on Nuclear Science, 21(3):21-42, June 1974.
[doi](https://doi.org/10.1109/TNS.1974.6499235)

This (default) version has low soft tissue contrast,
so it usually should be displayed
with a narrow window the using `clim` option:
=#

image1 = shepp_logan(256) # CT version by default
jim(image1, "SheppLogan (original CT version)", clim=(0.95, 1.05))


#=
## Over-sampling

When generating ellipse phantoms,
it is generally preferable
to "over-sample" the ellipse values within each pixel
to account for partial volume effect.
A factor of 3 over-sampling (along both axes) typically suffices,
so this factor is the default for the `shepp_logan` method.
Here is how the phantom image looks without over-sampling:
=#

image1o = shepp_logan(256; oversample=1)
jim(image1o, "No over sampling", clim=(0.95,1.05))


#=
Note that boundaries of the interior ellipses look smoother
in the original version with the default 3× over-sampling.
Most of the remaining examples use the recommended default over-sampling.
=#


#=
## Units

The original Shepp-Logan phantom
in the paper as cited above
did not give any units
for the grayscale values,
but the values were described
as "densities" and indeed the values
are reasonable for head CT scans
if interpreted as having g/cm³ units.
Here is a version with units.
=#

image1u = image1 * g/cm^3
jim(image1u, "SheppLogan with density units", clim=(0.95,1.05))


#=
## Toft version

This version is from:
Toft, Peter Aundal & Sørensen, John Aasted
"The Radon transform-theory and implementation,"
Technical University of Denmark (DTU), 1996. Page 201.
[pdf](https://files.openpdfs.org/1ra51GP6gJO.pdf)
=#

image2 = shepp_logan(256, SheppLoganToft())
jim(image2, "SheppLoganToft")


#=
## Emission tomography version

This version has low intensity for the skull
because typical PET/SPECT radiotracers
do not accumulate in bone regions.
It is probably also useful for MRI,
because typical MRI scans
have low signal from bone.
=#

image3 = shepp_logan(256, SheppLoganEmis())
jim(image3, "SheppLoganEmis")


#=
## BrainWeb version

This version was inspired by the
[BrainWeb](https://brainweb.bic.mni.mcgill.ca)
phantoms that have integer indices
for each of the different regions.
It should not be used directly,
but rather one should assign meaningful intensity values
to each of the integer indices.
=#

image4 = shepp_logan(256, SheppLoganBrainWeb())
jim(image4, "SheppLoganBrainWeb")


# For the BrainWeb version, there is no over-sampling by default,
# to preserve the integer indices.


#=
## Disjoint middle ellipses

Sometimes it can be more convenient
to have the middle ellipses be non-overlapping:
=#

params = ellipse_parameters(SheppLoganBrainWeb(), disjoint=true)
params = [(p[1:5]..., i) for (i, p) in enumerate(params)]
ob = ellipse(params)
x = range(-0.4, 0.4, 206)
y = range(-0.5, 0.5, 256)
oversample = 3
image5 = phantom(x, y, ob, oversample)
jim(x, y, image5, "Disjoint"; aspect_ratio = 1)


# ## Comedy version

image6 = shepp_logan(256, SouthPark(); fovs=(1,1))
jim(image6, "SouthPark")
