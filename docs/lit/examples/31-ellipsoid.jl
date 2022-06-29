#---------------------------------------------------------
# # [Ellipsoid](@id 31-ellipsoid)
#---------------------------------------------------------

#=
This page illustrates the `Ellipse` shape in the Julia package
[`ImagePhantoms`](https://github.com/JuliaImageRecon/ImagePhantoms.jl).

This page was generated from a single Julia file:
[31-ellipsoid.jl](@__REPO_ROOT_URL__/31-ellipsoid.jl).
=#

#md # In any such Julia documentation,
#md # you can access the source code
#md # using the "Edit on GitHub" link in the top right.

#md # The corresponding notebook can be viewed in
#md # [nbviewer](http://nbviewer.jupyter.org/) here:
#md # [`31-ellipsoid.ipynb`](@__NBVIEWER_ROOT_URL__/31-ellipsoid.ipynb),
#md # and opened in [binder](https://mybinder.org/) here:
#md # [`31-ellipsoid.ipynb`](@__BINDER_ROOT_URL__/31-ellipsoid.ipynb).


# ### Setup

# Packages needed here.

using ImagePhantoms: Ellipsoid, phantom, radon, spectrum
using ImageGeoms: ImageGeom, axesf
using MIRTjim: jim, prompt, mid3
using FFTW: fft, fftshift
using Unitful: mm, unit, °
using UnitfulRecipes
using Plots: plot, plot!, scatter!, gui, default
default(markerstrokecolor=:auto)


# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);

# ### Overview

#=
A basic shape used in constructing 3D digital image phantoms
is the ellipsoid, specified by its center, radii, angle(s) and value.
All of the methods in `ImagePhantoms` support physical units,
so we use such units throughout this example.
(Using units is recommended but not required.)

Here are 3 ways to define an ellipsoid object,
using parameters specified with physical units.
=#

if false # todo

center = (20mm, 10mm, 5mm)
radii = (35mm, 25mm, 15mm)
angles = (π/6, 0)
ob = Ellipsoid([40mm, 20mm, 10mm, 35mm, 25mm, 15mm, π/6, 0, 1.0f0])
ob = Ellipsoid(40mm, 20mm, 10mm, 35mm, 25mm, 15mm, π/6, 0, 1.0f0)
ob = Ellipsoid(center, radii, angles, 1.0f0) # recommended use


#=
### Phantom image using `phantom`

Make a 3D digital image of it using `phantom` and display it.
We use `ImageGeoms` to simplify the indexing.
=#

deltas = (1.2mm, 1.1mm, 1mm)
dims = (2^8, 2^8+2, 2^6)
offsets = (0.5, 0.5, 0.5) # for FFT spectra later
#x = (-M÷2:M÷2-1) * dx
#y = (-N÷2:N÷2-1) * dy
#jim(x, y, img)

ig = ImageGeom( ; dims, deltas, offsets)
img = phantom(axes(ig)..., [ob])
p1 = jim(axes(ig), img, "Ellipsoid phantom", xlabel="x", ylabel="y")


# Show middle slices

jim(mid3(img), "Middle 3 planes")

#=
### Spectrum using `spectrum`

There are two ways to examine the spectrum of this 3D image:
* using the analytical Fourier transform of the ellipse via `spectrum`
* applying the DFT via FFT to the digital image.
Because the shape has units `mm`, the spectra axes have units cycles/mm.
Appropriate frequency axes for DFT are provided by `axesf(ig)`.
=#

zscale = 1 / 4π / prod(radii) # normalize spectra by volume
spectrum_exact = spectrum(axesf(ig)..., [ob]) * zscale
sp = z -> max(log10(abs(z)/oneunit(abs(z))), -6) # log-scale for display
clim = (-6, 0) # colorbar limit for display
(xlabel, ylabel) = ("ν₁", "ν₂")
p2 = jim(axesf(ig), sp.(spectrum_exact), "log10|Spectrum|"; clim, xlabel, ylabel)


# Sadly `fft` cannot handle units currently, so this function is a work-around:
function myfft(x)
    u = unit(eltype(x))
    return fftshift(fft(fftshift(x) / u)) * u
end

spectrum_fft = myfft(img) * prod(ig.deltas) * zscale
p3 = jim(axesf(ig), sp.(spectrum_fft), "log10|DFT|"; clim, xlabel, ylabel)


# Compare the DFT and analytical spectra to validate the code
@assert maximum(abs, spectrum_exact - spectrum_fft) /
        maximum(abs, spectrum_exact) < 3e-3
p4 = jim(axesf(ig), abs.(spectrum_fft - spectrum_exact), "Difference"; xlabel, ylabel)
jim(p1, p4, p2, p3)


#=
### Parallel-beam projections using `radon`

Compute 2D projection views of the object using `radon`.
Validate it using the projection-slice theorem aka Fourier-slice theorem.
=#

end # todo

du, dv = 0.3mm, 0.5mm # projection view pixel spacing
nu, nv = 2^9, 2^8 # projection view dimensions
u = (-nu÷2:nu÷2-1) * du # projection view samples
v = (-nv÷2:nv÷2-1) * dv
fu = (-nu÷2:nu÷2-1) / nu / du # corresponding spectral axes
fv = (-nv÷2:nv÷2-1) / nv / dv
ϕs, θs = :(π/2), :(π/7)
ϕ, θ = eval.((ϕs, θs))
proj = radon(u, v, [ϕ], [θ], [ob]) # single projection of a single object
#p5 = jim(u, v, proj, "Projection at ϕ=$ϕs θ=$θs", xlabel="u", ylabel="v")


ϕd = 0:6:360
ϕs = deg2rad.(ϕd) # * Unitful.rad # todo round unitful Unitful.°
proj = radon(u, v, ϕs, [θ], [ob]) # many projection views
#src p5 = jim(u, v, proj; title="projection views")

clim = 2 * maximum(radii) * ob.value; clim = (zero(clim), clim)
for (ip, ϕ) in enumerate(ϕd)
    jim(u, v, proj[:,:,ip,1]; title="ϕ=$ϕ", clim, prompt=false)
    gui() # todo anim
end


#=
The above sampling generated a parallel-beam projection,
but one could make a cone-beam projection
by sampling `(u, v, ϕ, θ)` appropriately.
See Sinograms.jl.
=#

todo


# ### Fourier-slice theorem illustration

# Pick one particular view angle (55°) and look at its slice and spectra.
ia = argmin(abs.(ϕ .- 55°))
slice = sino[:,ia]
slice_fft = myfft(slice) * dr
angle = round(rad2deg(ϕ[ia]), digits=1)

kx, ky = (fr * cos(ϕ[ia]), fr * sin(ϕ[ia])) # Fourier-slice theorem
slice_ft = spectrum(ob).(kx, ky)
@assert maximum(abs, slice_ft - slice_fft) / maximum(abs, slice_ft) < 2e-4

p3 = plot(r, slice, title="profile at ϕ = $angle", label="")
p4 = plot(title="1D spectra")
scatter!(fr, abs.(slice_fft), label="abs fft", color=:blue)
scatter!(fr, real(slice_fft), label="real fft", color=:green)
scatter!(fr, imag(slice_fft), label="imag fft", color=:red, xlims=(-1,1).*(1.0/mm))

plot!(fr, abs.(slice_ft), label="abs", color=:blue)
plot!(fr, real(slice_ft), label="real", color=:green)
plot!(fr, imag(slice_ft), label="imag", color=:red)
plot(p1, p5, p3, p4)


#=
The good agreement between the analytical spectra (solid lines)
and the DFT samples (disks)
validates that `phantom`, `radon`, and `spectrum`
are all self consistent for this `Ellipsoid` object.
=#



#=
todo: pieces below here
=#

    # test sinogram with projection-slice theorem

    du,dv = 0.02m, 0.03m
    nu,nv = 2^9, 2^8
    u = (-nu÷2:nu÷2-1) * du
    v = (-nv÷2:nv÷2-1) * dv
    fu = (-nu÷2:nu÷2-1) / nu / du
    fv = (-nv÷2:nv÷2-1) / nv / dv
    ϕ = deg2rad.(0:6:180) # * Unitful.rad # todo round unitful Unitful.°
    θ = [π/7]
    sino = @inferred radon(u, v, ϕ, θ, [ob])

#=
todo
    ia = argmin(abs.(ϕ .- deg2rad(55)))
    slice = sino[:,ia]
    Slice = myfft(slice) * dr
    angle = round(rad2deg(ϕ[ia]), digits=1)

    kx, ky = (fr * cos(ϕ[ia]), fr * sin(ϕ[ia])) # Fourier-slice theorem
    ideal = spectrum(ob).(kx, ky)

if DEBUG
    p2 = jim(r, rad2deg.(ϕ), sino; aspect_ratio=:none, title="sinogram")
    jim(p1, p2)
    p3 = plot(r, slice, title="profile at ϕ = $angle", label="")
    p4 = scatter(fr, abs.(Slice), label="abs fft", color=:blue)
    scatter!(fr, real(Slice), label="real fft", color=:green)
    scatter!(fr, imag(Slice), label="imag fft", color=:red,
        xlims=(-1,1).*(1.2/m), title="1D spectra")

    plot!(fr, abs.(ideal), label="abs", color=:blue)
    plot!(fr, real(ideal), label="real", color=:green)
    plot!(fr, imag(ideal), label="imag", color=:red)
    plot(p1, p2, p3, p4); gui()
end

    @test maximum(abs, ideal - Slice) / maximum(abs, ideal) < 2e-4
=#
