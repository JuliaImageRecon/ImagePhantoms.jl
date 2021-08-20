#---------------------------------------------------------
# # [Shepp-Logan Phantoms](@id 8-shepp)
#---------------------------------------------------------

# This page illustrates the Shepp-Logan phantoms in the Julia package
# [`ImagePhantoms`](https://github.com/JuliaImageRecon/ImagePhantoms.jl).

# ### Setup

# Packages needed here.

using ImagePhantoms
using ImageGeoms: ImageGeom, axesf
using MIRTjim: jim, prompt

# The following line is helpful when running this example.jl file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);

# ### Overview

# There are several variations of the Shepp-Logan phantom available.


# ### CT version

# Original version from:
# Larry A Shepp, Benjamin F Logan,
# "The Fourier reconstruction of a head section,"
# IEEE Transactions on Nuclear Science, 21(3):21-42, June 1974.
# [doi](http://doi.org/10.1109/TNS.1974.6499235)
#
# This (default) version has low soft tissue contrast,
# so it usually should be displayed
# with a narrow window the using `clim` option:

image1 = shepp_logan(256) # CT version by default
jim(image1, "SheppLogan (original CT version)", clim=(0.95, 1.05), yflip=false)


# ### Over-sampling

# When generating ellipse phantoms,
# it is generally preferable
# to "over-sample" the ellipse values within each pixel
# to account for partial volume effect.
# A factor of 3 over-sampling (along both axes) typically suffices,
# so this factor is the default for the `shepp_logan` method.
# Here is how the phantom image looks without over-sampling:

image1o = shepp_logan(256; oversample=1)
jim(image1o, "No over sampling", clim=(0.95,1.05), yflip=false)


# Notice that boundaries of the interior ellipses look smoother
# in the original version with the default 3× over-sampling.
# Most of the remaining examples use the recommended default over-sampling.


# ### Toft version

# This version is from:
# Toft, Peter Aundal & Sørensen, John Aasted
# "The Radon transform-theory and implementation,"
# Technical University of Denmark (DTU), 1996. Page 201.
# [pdf](https://files.openpdfs.org/1ra51GP6gJO.pdf)

image2 = shepp_logan(256, SheppLoganToft())
jim(image2, "SheppLoganToft", yflip=false)


# ### Emission tomography version

# This version has low intensity for the skull
# because typical PET/SPECT radiotracers
# do not accumulate in bone regions.
# It is probably also useful for MRI,
# because typical MRI scans
# have low signal from bone.

image3 = shepp_logan(256, SheppLoganEmis())
jim(image3, "SheppLoganEmis", yflip=false)


# ### BrainWeb version

# This version was inspired by the
# [BrainWeb](https://brainweb.bic.mni.mcgill.ca)
# phantoms that have integer indices
# for each of the different regions.
# It should not be used directly,
# but rather one should assign meaningful intensity values
# to each of the integer indices.

image4 = shepp_logan(256, SheppLoganBrainWeb())
jim(image4, "SheppLoganBrainWeb", yflip=false)


# For the BrainWeb version, there is no over-sampling by default,
# to preserve the integer indices.


# ### Comedy version

image5 = shepp_logan(256, SouthPark(); fovs=(1,1), oversample=3)
jim(image5, "SouthPark", yflip=false)
