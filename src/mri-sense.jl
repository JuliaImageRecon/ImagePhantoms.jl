#=
mri-sense.jl
Methods for generating phantom spectra based on the 2012 paper
[Guerquin-Kern et al.](http://doi.org/10.1109/TMI.2011.2174158)
that combines analytical k-space values of the phantom
with an analytical model for the sensitivity maps.
=#

using LazyGrids: ndgrid
using ImagePhantoms: spectrum, Object2d

export smap_basis, smap_fit, spectrum, spectra


"""
   smap_basis(mask ; kmax, kt, ki)

Construct Fourier basis for representing MRI sensitivity maps
in terms of separable complex exponential signals, products of
of the form `basis = (k,N) -> exp.(2im * π * (0:N-1) * kfun(k,N))`.
The default is `kfun(k,N) = k / 2N`, which is DCT-II like frequencies,
leading to better boundary behavior than the DFT frequencies `k/N`.

# Input
* `mask::AbstractArray{Bool,D}` binary support mask for region to reconstruct

# Option
* `kmax::Int = 5` default frequency index -kmax:kmax in all dimensions
* `kfun::Function = (k,N) -> k / (2N)` # DCT-II frequency
* `deltas::NTuple{D,<:Number} = ones(D)` pixel sizes
(For additional options `kmaxs`, `kt`, `ki`, see code.)

# Output
* `(; B, ν)` where `B` is basis matrix of size `count(mask) × nk`
  where typically `nk = (2*kmax+1)^D`
  and `ν` is `nk` frequency tuples;
  each tuple has form `ν = kfun.(Tuple(k), size(mask)) ./ deltas`.
"""
function smap_basis(
    mask::AbstractArray{Bool,D} ;
    kmax::Int = 5,
    kmaxs::NTuple{D, Int} = ntuple(i -> kmax, D),
    kt::NTuple{D, UnitRange{<:Integer}} = ntuple(i -> -kmaxs[i]:kmaxs[i], D),
    ki::CartesianIndices{D} = CartesianIndices(kt),
    deltas::NTuple{D,<:Number} = ones(D),
    kfun::Function = (k,N) -> k / (2N), # DCT-II frequency
    basis::Function = (k,N) -> exp.(2im * π * (0:N-1) * kfun(k,N)), # 1D basis
    T::DataType = ComplexF32,
) where D
@show kmax
    B = zeros(T, count(mask), length(ki))
    for (i, k) in enumerate(ki)
        tmp = ndgrid([basis(k[d], size(mask,d)) for d = 1:D]...)
        tmp = .*(tmp...)
        B[:,i] = tmp[mask]
    end

    ν = [kfun.(Tuple(k), size(mask)) ./ deltas for k in ki]
    return (; B, ν)
end


"""
    smap_fit(smaps, embed ; mask, kwargs...)

Fit MRI sensitivity maps `smaps` using `smap_basis(mask ; kwargs...)`.
Provide `ImageGeoms.embed` or equivalent.
Return named tuple `(B, ν, coefs, nrmse, smaps)`:
* `(B, ν)` from `smap_basis`
* `coefs::Vector` : `[ncoil]` each of length `nk`
* `nrmse::Real` : `smaps` vs `smaps_fit`
* `smaps::Vector{Array{D}}` : `smaps_fit`
"""
function smap_fit(
    smaps::Vector{<:AbstractArray{<:Number,D}},
	embed::Function ;
    mask::AbstractArray{Bool,D} = trues(size(smaps[1])),
    kwargs...
) where D
    (B, ν) = smap_basis(mask ; kwargs...)
    coefs = [B \ smap[mask] for smap in smaps] # coefficients for each smap

    smaps_fit = [embed(B * coef, mask) for coef in coefs]
    ss(v) = sum(x -> sum(abs2, x), v)
    @show nrmse = sqrt( ss(smaps_fit - smaps) / ss(smaps) )
    return (; B, ν, coefs, nrmse, smaps = smaps_fit)
end


"""
    spectrum(ob::Object2d, coefs::AbstractVector, f::)
Version of `spectrum(ob)` suitable for parallel MRI with sensitivity maps
that were fit previously using `smap_fit` for a single coil
with fit coefficients `coefs` and frequencies `f` (vector of tuples).
"""
function spectrum(ob::Object2d, coefs::AbstractVector, f::Any)
    return (fx,fy) -> sum(
        i -> coefs[i] * spectrum(ob)(fx - f[i][1], fy - f[i][2]),
        1:length(f),
    )
end


"""
    spectrum(ob::Object2d, fit::NamedTuple, coil::Int)
Version of `spectrum(ob)` suitable for parallel MRI with sensitivity maps
that were fit previously using `smap_fit`.
"""
spectrum(ob::Object2d, fit::NamedTuple, coil::Int) =
    spectrum(ob, fit.coefs[coil], fit.ν)


"""
    spectra(fx, fy, oa::Array{<:Object2d}, fit::NamedTuple)
Version of `spectrum` suitable for parallel MRI with sensitivity maps
that were fit previously using `smap_fit`.
Returns a Vector of `ncoil` kspace data Vectors of dimension `length(fx)`
"""
function spectra(
    fx::AbstractVector,
    fy::AbstractVector,
    oa::Array{<:Object2d},
    fit::NamedTuple,
)
    ncoil = length(fit.smaps)
    return [sum(ob -> spectrum(ob, fit, ic).(fx,fy), oa) for ic=1:ncoil]
end
