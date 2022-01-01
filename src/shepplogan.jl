#=
shepplogan.jl
=#

export shepp_logan, ellipse_parameters
export SheppLogan
export SheppLoganEmis
export SheppLoganBrainWeb
export SheppLoganToft
export SouthPark


"""
    EllipsePhantomVersion
Parent type for different versions of ellipse phantoms:
* `SheppLogan` original CT version from Shepp&Logan paper
* `SheppLoganEmis` higher contrast version suitable for emission tomography
* `SheppLoganBrainWeb` integer index version based on brainweb
* `SheppLoganToft` higher contrast version from Toft, 1996
* `SouthPark` for fun
"""
abstract type EllipsePhantomVersion end

"""
    SheppLogan
Original version from :
Larry A Shepp, Benjamin F Logan,
"The Fourier reconstruction of a head section,"
IEEE Transactions on Nuclear Science, 21(3):21-42, June 1974.
[doi](http://doi.org/10.1109/TNS.1974.6499235)

Also in Kak and Slaney 1988 text, p. 255.
[doi](http://doi.org/10.1137/1.9780898719277)
"""
struct SheppLogan <: EllipsePhantomVersion end

"""
    SheppLoganToft
Toft, Peter Aundal & Sørensen, John Aasted
"The Radon transform-theory and implementation,"
Technical University of Denmark (DTU), 1996. Page 201.
https://files.openpdfs.org/1ra51GP6gJO.pdf
"""
struct SheppLoganToft <: EllipsePhantomVersion end

struct SheppLoganEmis <: EllipsePhantomVersion end
struct SheppLoganBrainWeb <: EllipsePhantomVersion end
struct SouthPark <: EllipsePhantomVersion end


"""
    values = shepp_logan_values(::EllipsePhantomVersion)
Return 10 Shepp-Logan ellipse amplitudes for various versions.
"""
shepp_logan_values(::SheppLogan) =
    [2.0, -0.98, -0.02, -0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
shepp_logan_values(::SheppLoganEmis) =
    [1, 1, -2, 2, 3, 4, 5, 6, 1, 1]
shepp_logan_values(::SheppLoganBrainWeb) =
    [1, 0, 2, 3, 4, 5, 6, 7, 8, 9] # brainweb uses index 1-10
shepp_logan_values(::SheppLoganToft) =
    [1.0, -0.8, -0.2, -0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]


"""
    phantom = Ellipse(n × 6 AbstractMatrix)
Return vector of Ellipse objects, one for each row of input matrix.
"""
function Ellipse(params::AbstractMatrix{<:RealU})
    size(params,2) == 6 || throw("ellipses need 6 parameters")
    return [Ellipse(params[n,:]) for n in 1:size(params,1)]
end

function shepp_logan(case::EllipsePhantomVersion; kwargs...)
    return Ellipse(ellipse_parameters(case; kwargs...))
end


"""
    image = shepp_logan(M, [N,], case, options...)
Convenience method for generating `M×N` samples of Shepp-Logan phantoms.

In
* `M::Int` : horizontal size
* `N::Int` : vertical size, defaults to `M`
* `case::EllipsePhantomVersion = SheppLogan()`

# Options
* `oversample::Int = 3` (usually)
* `yflip::Bool = true` (reverse `y` samples for convenience.)
* `kwargs...` remaining options passed to `ellipse_parameters` for parameters.

# Out
* `image` : `M × N` matrix

The default here is 3× over-sampling along both axes (9 samples per pixel),
except for the `SheppLoganBrainWeb` phantom that consists of integer indices.
"""
function shepp_logan(
    M::Int, N::Int,
    case::EllipsePhantomVersion = SheppLogan() ;
    oversample::Int = case == SheppLoganBrainWeb() ? 1 : 3,
    yflip::Bool = true,
    kwargs...
)
    ob = shepp_logan(case ; kwargs...)
    x = LinRange(-0.5, 0.5, M)
    y = LinRange(-0.5, 0.5, N)
    if yflip
        y = reverse(y)
    end
    return oversample > 1 ?
        phantom(x, y, ob, oversample) :
        phantom(x, y, ob) # type unstable because of over-sampling
end

shepp_logan(M::Int, case::EllipsePhantomVersion = SheppLogan(); kwargs...) =
    shepp_logan(M, M, case ; kwargs...)


"""
    params = ellipse_parameters(case::Symbol; fovs::NTuple{2}, u::Tuple, disjoint::Bool)

By default the first four columns are unitless "fractions of field of view",
so columns 1,3 are scaled by `xfov` and columns 2,4 are scaled by `yfov`,
where `(xfov, yfov) = fovs`.
The optional 3-tuple `u` specifies scaling and/or units:
* columns 1-4 (center, radii) are scaled by `u[1]` (e.g., mm),
* column 5 (angle) is scaled by `u[2]` (e.g., `1` or `°`),
* column 6 (value) is scaled by `u[3]` (e.g., `1/cm`) for an attenuation map.
If `disjoint==true` then the middle ellipse positions are adjusted to avoid overlap.
"""
function ellipse_parameters(
    case::EllipsePhantomVersion = SheppLogan() ;
    fovs::NTuple{2,RealU} = (1,1),
    u::NTuple{3,Any} = (1,1,1), # unit scaling
    disjoint::Bool = false,
)
    (xfov, yfov) = fovs

    # The "Number" here is to enable units
    params = Number[ # original CT version
    0       0       0.92    0.69    90    2     # skull
    0       -0.0184 0.874   0.6624  90    -0.98 # brain
    0.22    0       0.31    0.11    72    -0.02 # right big
    -0.22   0       0.41    0.16   108    -0.02 # left big
    0       0.35    0.25    0.21    90    0.01 # top
    0       0.1     0.046   0.046    0    0.01 # middle high
    0       -0.1    0.046   0.046    0    0.01 # middle low
    -0.08   -0.605  0.046   0.023    0    0.01 # bottom left
    0       -0.605  0.023   0.023    0    0.01 # bottom center
    0.06    -0.605  0.046   0.023   90    0.01 # bottom right
    ]

    if disjoint # tweak so that middle ellipses are disjoint
        params[4,1] += -0.04
        params[5,2] += 0.07
    end

    params[:,[1,3]] .*= xfov/2
    params[:,[2,4]] .*= yfov/2
    params[:,5] .*= π/180
    params[:,6] = shepp_logan_values(case)

    return Number[params[:,1:4] * u[1] params[:,5] * u[2] params[:,6] * u[3]]
end


"""
    params = ellipse_parameters(SouthPark() ; fovs::NTuple{2,Number} = (100,100))
Ellipse parameters for "South Park" phantom.
"""
function ellipse_parameters(::SouthPark ; fovs::NTuple{2,RealU} = (100,100))
    params = Float64[
        0 0 85 115 0 100
        0 -60 30 20 0 -80 # mouth
        30 20 25 35 30 20 # eyes
        -30 20 25 35 -30 20
        35 25 7 7 0 -100 # pupils
        -15 25 7 7 0 -100
        0 75 60 15 0 -50  # hat
    ]
    params[:,[1,3]] .*= fovs[1] / 256
    params[:,[2,4]] .*= fovs[2] / 256
    params[:,5] .*= π/180
    return params
end
