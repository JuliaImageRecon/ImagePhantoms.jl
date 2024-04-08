#=
shepplogan.jl
=#

export shepp_logan, ellipse_parameters
export ellipsoid_parameters
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
[doi](https://doi.org/10.1109/TNS.1974.6499235)

Also in Kak and Slaney 1988 text, p. 255.
[doi](https://doi.org/10.1137/1.9780898719277)
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
    phantom = ellipse(Vector{Tuple{6 params}})
Return vector of `Object{Ellipse}`,
one for each element of input vector of tuples.
Often the input comes from `ellipse_parameters`.
"""
function ellipse(params::Vector{<:Tuple})
    length(params[1]) == 6 || throw("ellipses need 6 parameters")
    out = [ellipse(p...) for p in params]
    return out
end


function shepp_logan(case::EllipsePhantomVersion; kwargs...)
    return ellipse(ellipse_parameters(case; kwargs...))
end


"""
    image = shepp_logan(M, [N,] [case] ; options...)
Convenience method for generating `M×N` samples of Shepp-Logan phantoms.

# In
* `M::Int` : horizontal size
* `N::Int` : vertical size, defaults to `M`
* `case::EllipsePhantomVersion = SheppLogan()`

# Options
* `oversample::Int = 3` (usually)
* `yflip::Bool = true` (reverse `y` samples for convenience.)
* `T::Type{<:Number}` default `Float32` (except `Int` for `BrainWeb` version)
* `kwargs...` remaining options passed to `ellipse_parameters` for parameters.

# Out
* `image` : `M × N` matrix

The default here is 3× over-sampling along both axes (9 samples per pixel),
except for the `SheppLoganBrainWeb` phantom that consists of integer indices.
"""
function shepp_logan(
    M::Int, N::Int,
    case::EllipsePhantomVersion ;
    oversample::Int = 3,
    yflip::Bool = true,
    T::Type{<:Number} = Float32,
    kwargs...
)
    ob = shepp_logan(case ; kwargs...)
    x = range(-0.5, 0.5, M)
    y = range(-0.5, 0.5, N)
    if yflip
        y = reverse(y)
    end
    out = oversample > 1 ?
        phantom(x, y, ob, oversample ; T) :
        T.(phantom(x, y, ob))
    return out::Matrix{T}
end


# For type stability, BrainWeb case (with integer values) is separate
function shepp_logan(
    M::Int, N::Int,
    case::SheppLoganBrainWeb ;
    oversample::Int = 1,
    yflip::Bool = true,
    kwargs...
)
    oversample == 1 || throw("oversample must be 1 for BrainWeb()")
    ob = shepp_logan(case ; kwargs...)
    x = range(-0.5, 0.5, M)
    y = range(-0.5, 0.5, N)
    if yflip
        y = reverse(y)
    end
    out = Int.(phantom(x, y, ob))
    return out
end


shepp_logan(M::Int, case::EllipsePhantomVersion = SheppLogan(), args...; kwargs...) =
    shepp_logan(M, M, case, args... ; kwargs...)

shepp_logan(M::Int, N::Int ; kwargs...) =
    shepp_logan(M, N, SheppLogan() ; kwargs...)


"""
    ellipse_parameters_shepplogan( ; disjoint::Bool)
`10 × 6 Matrix{Float64}` of classic Shepp-Logan ellipse parameters.
If `disjoint==true` then the middle ellipse positions are adjusted to avoid overlap.
"""
function ellipse_parameters_shepplogan( ; disjoint::Bool = false)

    params = Float64[ # original CT version
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

    return params
end


"""
    ellipse_parameters_uscale(params, fovs, uc, ua, uv)
Return vector of Tuples after FOV and unit scaling.
"""
function ellipse_parameters_uscale(
    params::Matrix{Tp},
    fovs::NTuple{2,Tf},
    uc::Tc,
    ua::Ta,
    uv::Tv,
) where {
    Tp <: AbstractFloat,
    Tf <: RealU,
    Tc <: RealU,
    Ta <: RealU,
    Tv <: Number,
}

    C = typeof(oneunit(Tf) * oneunit(Tc) * one(Tp))
    A = typeof(oneunit(Ta) * one(Tp))
    V = typeof(oneunit(Tv) * one(Tp))
    T = Tuple{C, C, C, C, A, V}

    size(params,2) == 6 || throw(ArgumentError("params not N × 6"))
    N = size(params,1)
    out = Vector{T}(undef, N)
    for n in 1:N
        tmp = params[n,:]
        tmp = (
            (tmp[1:2] * uc .* fovs)...,
            (tmp[3:4] * uc .* fovs)...,
            tmp[5] * ua,
            tmp[6] * uv,
        )
        out[n] = tmp
    end
    return out
end


"""
    ellipse_parameters(case; fovs::NTuple{2}, u::NTuple{3}, disjoint::Bool)

Return vector of Tuples of ellipse parameters.
By default the first four elements of each tuple
are unitless "fractions of field of view",
so elements 1,3 are scaled by `xfov`
and elements 2,4 are scaled by `yfov`,
where `(xfov, yfov) = fovs`.
The optional 3-tuple `u` specifies scaling and/or units:
* elements 1-4 (center, radii) are scaled by `u[1]` (e.g., mm),
* elements 5 (angle) is scaled by `u[2]` (e.g., `1` or `°`),
* elements 6 (value) is scaled by `u[3]` (e.g., `1/cm`) for an attenuation map.
If `disjoint==true` then the middle ellipse positions are adjusted to avoid overlap.
"""
function ellipse_parameters(
    case::EllipsePhantomVersion = SheppLogan() ;
    fovs::NTuple{2,RealU} = (1,1),
    u::NTuple{3,Number} = (1,1,1), # unit scaling
    disjoint::Bool = false,
)

    params = ellipse_parameters_shepplogan( ; disjoint)
    params[:,1:4] ./= 2
    params[:,5] .*= π/180
    params[:,6] = shepp_logan_values(case)

    out = ellipse_parameters_uscale(params, fovs, u[1], u[2], u[3])
    return out
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

    params[:,[1,3]] ./= 256
    params[:,[2,4]] ./= 256
    params[:,5] .*= π/180

    out = ellipse_parameters_uscale(params, fovs, 1, 1, 1)
    return out
end


# 3D

abstract type EllipsoidPhantomVersion end

struct SheppLogan3 <: EllipsoidPhantomVersion end


"""
    ellipsoid_parameters_shepplogan( ; disjoint::Bool)
`12 × 10 Matrix{Float64}` of 3D Shepp-Logan ellipsoid parameters
`(cx,cy,cz, rx,ry rz, Φ,Θ=0,Ψ=0, ρ)`.
By default the first 6 columns are unitless "fractions of field of view".
"""
function ellipsoid_parameters_shepplogan()

#=
The following 3D ellipsoid parameters came from leizhu@stanford.edu
who said that the Kak&Slaney 1988 values are incorrect.
=#
#       x       y       z       rx      ry      rz      Φ°     density
    params = Float64[
        0       0       0       0.69    0.92    0.9     0       2.0;
        0       -0.0184 0       0.6624  0.874   0.88    0       -0.98;
        -0.22   0       -0.25   0.41    0.16    0.21    -72     -0.02;
        0.22    0       -0.25   0.31    0.11    0.22    72      -0.02;
        0       0.35    -0.25   0.21    0.25    0.35    0       0.01;
        0       0.1     -0.25   0.046   0.046   0.046   0       0.01;
        -0.08   -0.605  -0.25   0.046   0.023   0.02    0       0.01;
        0       -0.1    -0.25   0.046   0.046   0.046   0       0.01;
        0       -0.605  -0.25   0.023   0.023   0.023   0       0.01;
        0.06    -0.605  -0.25   0.046   0.023   0.02    -90     0.01;
        0.06    -0.105  0.0625  0.056   0.04    0.1     -90     0.02;
        0       0.1     0.625   0.056   0.056   0.1     0       -0.02
    ]
    params[:,1:6] ./= 2 # radii
    params[:,7] .*= π/180 # radians
    Θ = zeros(size(params,1))
    Ψ = zeros(size(params,1))
    out = [params[:,1:7] Θ Ψ params[:,8]] # (12,10)
    return out
end


"""
    oa = ellipsoid(Vector{Tuple{10 params}})
Return vector of `Object{Ellipsoid}`,
one for each element of input vector of tuples.
Often the input comes from `ellipsoid_parameters`.
"""
function ellipsoid(params::Vector{<:Tuple})
    length(params[1]) == 10 || throw("ellipsoids need 10 parameters")
    out = [ellipsoid(p...) for p in params]
    return out
end


"""
    ellipsoid_parameters_uscale(params, fovs, uc, ua, uv)
Return vector of Tuples after FOV and unit scaling.
"""
function ellipsoid_parameters_uscale(
    params::Matrix{Tp},
    fovs::NTuple{3,Tf},
    uc::Tc,
    ua::Ta,
    uv::Tv,
) where {
    Tp <: AbstractFloat,
    Tf <: RealU,
    Tc <: RealU,
    Ta <: RealU,
    Tv <: Number,
}
    C = typeof(oneunit(Tf) * oneunit(Tc) * one(Tp))
    A = typeof(oneunit(Ta) * one(Tp))
    V = typeof(oneunit(Tv) * one(Tp))
    T = Tuple{C, C, C, C, C, C, A, A, A, V}

    size(params,2) == 10 || throw(ArgumentError("params not N × 10"))
    N = size(params,1)
    out = Vector{T}(undef, N)
    for n in 1:N
        tmp = params[n,:]
        tmp = (
            (tmp[1:3] * uc .* fovs)...,
            (tmp[4:6] * uc .* fovs)...,
            (tmp[7:9] * ua)...,
            tmp[10] * uv,
        )
        out[n] = tmp
    end
    return out
end


"""
    ellipsoid_parameters(case; fovs::NTuple{3}, u::NTuple{3})

Return vector of Tuples of ellipsoid parameters.
By default the parameters are for a 3D Shepp-Logan ellipsoid phantom.
By default the first 6 elements of each tuple
are unitless "fractions of field of view",
so elements 1,4 are scaled by `xfov`
and elements 2,5 are scaled by `yfov`,
and elements 3,6 are scaled by `zfov`,
where `(xfov, yfov, zfov) = fovs`.
The optional 3-tuple `u` specifies scaling and/or units:
* elements 1-6 (center, radii) are scaled by `u[1]` (e.g., mm),
* elements 7-8 (angles) are scaled by `u[2]` (e.g., `1` or `°`),
* element 9 (value) is scaled by `u[3]` (e.g., `1/cm`) for an attenuation map.
"""
function ellipsoid_parameters(
    case::EllipsoidPhantomVersion = SheppLogan3() ;
    fovs::NTuple{3,RealU} = (1,1,1),
    u::NTuple{3,Number} = (1,1,1), # unit scaling
)

    params = ellipsoid_parameters_shepplogan() # (N,10)

    case == SheppLogan3() || error("unsupported case $case")
#   params[:,10] = shepp_logan_values(case)

    out = ellipsoid_parameters_uscale(params, fovs, u...)
end
