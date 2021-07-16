#=
shape.jl
Types for describing 2D and 3D shapes, like ellipses and ellipsoids.
=#

export AbstractShape, AbstractShape2, AbstractShape3


"""
    AbstractShape
Generic shape type for `ImagePhantoms`.
"""
abstract type AbstractShape end

"""
    AbstractShape2 <: AbstractShape
Generic 2D shape type (with subtypes Ellipse...)
"""
abstract type AbstractShape2 <: AbstractShape end

"""
    AbstractShape2 <: AbstractShape
Generic 3D shape type (with subtypes Ellipsoid...)
"""
abstract type AbstractShape3 <: AbstractShape end
