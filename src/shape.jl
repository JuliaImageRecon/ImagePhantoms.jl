#=
shape.jl
Types for describing 2D and 3D shapes, like ellipses and ellipsoids.
=#

export AbstractShape


"""
    AbstractShape{D}
Generic shape type for `ImagePhantoms`.
The dimension `D` is likely `2` or `3`,
e.g., `2` for an `Ellipse` and `3` for an `Ellipsoid`.
"""
abstract type AbstractShape{D} end
