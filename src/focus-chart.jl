#=
focus-chart.jl
generate "focus chart" phantom from triangles
2021-09-12, Jeff Fessler, University of Michigan
=#

export focus_chart

using ImagePhantoms: Triangle, triangle, Object2d, RealU


"""
    focus_chart( ; radius=1, nspoke=30, value=1)

Generate `nspoke` Triangle phantom parameters for a focus chart.
Returns `Vector{Object2d{Triangle}` for passing to `phantom`.

# Options
- `radius::RealU = 1` radius of phantom
- `nspoke::Int = 60` # of spokes
- `value::Number = 1` alternate between 0 and this value
"""
focus_chart( ;
    nspoke::Int = 60,
    radius::RealU = 1,
    value::Number = 1,
    values::AbstractVector = value * (iseven.(1:nspoke)*2 .- 1), # alternating
) = focus_chart.(1:nspoke, nspoke, radius, values)


# ith of nspoke triangles
function focus_chart(i::Int, nspoke::Int, radius::RealU, value::Number)
    ϕ = (i-1)/nspoke * 2π
    (s, c) = sincos(ϕ)
    θ = π/nspoke
    (w,h) = sincos(θ) .* radius .* (2, 2/sqrt(3))
    offset = cos(θ) * radius + eps() # trick to avoid center pile up
    return triangle((s, -c) .* offset, (w,h), ϕ, value)
end
