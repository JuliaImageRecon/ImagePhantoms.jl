#=
cylinder.jl
=#


export Cylinder, cylinder


"""
    Cylinder <: AbstractShape{3}
"""
struct Cylinder <: AbstractShape{3} end


# constructor


"""
    cylinder(cx, cy, cz, wx, wy, wz, Φ, Θ, value::Number)
    cylinder(center::NTuple{3,RealU}, width::NTuple{3,RealU}, angle::NTuple{3,RealU}, v)
Construct `Object{Cylinder}` from parameters;
here `width` is the *radius* in x,y and the *height* in z.
"""
cylinder(args... ; kwargs...) = Object(Cylinder(), args...; kwargs...)


# methods


volume1(::Cylinder) = π

ℓmax1(::Cylinder) = √5 # max line integral through unit cylinder

ℓmax(ob::Object3d{Cylinder}) = sqrt(sum(abs2, (2maximum(ob.width[1:2]), ob.width[3])))


"""
    phantom1(ob::Object3d{Cylinder}, (x,y,z))
Evaluate unit cylinder at `(x,y,z)`,
for unitless coordinates.
"""
phantom1(ob::Object3d{Cylinder}, xyz::NTuple{3,Real}) =
    (abs(xyz[3]) ≤ 0.5) && (sum(abs2, xyz[1:2]) ≤ 1)


# radon


# line integral through rectangle of width (wx,wy) at (r,ϕ)
function _rect_proj(wx::Real, wy::Real, r::Real, ϕ::Real)
    (sϕ, cϕ) = sincos(ϕ)
    rp = sqrt((wx * cϕ)^2 + (wy * sϕ)^2) # projected radius
    dis = r / rp # scaled distance from center
    abs_cos_ang_pi = wx * abs(cϕ) / rp
    abs_sin_ang_pi = wy * abs(sϕ) / rp

    # break points of the trapezoid
    len = 1 / max(abs_cos_ang_pi, abs_sin_ang_pi)
    dmax = (abs_cos_ang_pi + abs_sin_ang_pi) / 2
    dbreak = abs(abs_cos_ang_pi - abs_sin_ang_pi) / 2
    dmax_break = dmax - dbreak
    scale = wx * wy / rp * len #lmax

    return scale * trapezoid(dis, -dmax, -dbreak, dbreak, dmax)
end


# x-ray transform (line integral) of unit cylinder
# `u,v` should be unitless
function xray1(
    ::Cylinder,
    u::Real,
    v::Real,
    ϕ::RealU, # azim (irrelevant)
    θ::RealU, # polar
)
    T = promote_type(typeof(u), typeof(v), Float32)

    r = abs(u)
    if r > 1
        return zero(T)
    end

    # rectangle in plane of distance `r` from origin
    wz = 1 # from unit-height of cylinder
    wy = 2 * sqrt(1 - r^2)

    return T(_rect_proj(wz, wy, v, θ))
end


# spectrum

"""
    spectrum(::Object3d{Cylinder}, (kx,ky,kz))
Spectrum of unit cylinder at `(kx,ky,kz)`,
for unitless spatial frequency coordinates.
"""
function spectrum1(::Object3d{Cylinder}, kxyz::NTuple{3,Real})
    return 4 * jinc(2 * sqrt(sum(abs2, kxyz[1:2]))) * sinc(kxyz[3])
end
