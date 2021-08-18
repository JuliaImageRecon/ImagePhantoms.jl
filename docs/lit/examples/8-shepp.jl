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
jim(image1, "SheppLogan (original CT version)", clim=(0.9, 1.1), yflip=false)


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


# ### Comedy version

image5 = shepp_logan(256, SouthPark(); fovs=(1,1))
jim(image5, "SouthPark", yflip=false)
