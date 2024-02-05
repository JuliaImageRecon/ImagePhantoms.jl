#=
triangle.jl
2D triangle shape
=#

using ImagePhantoms #: Object2d

export Triangle

"""
    Triangle{T} <: AbstractShape2
By default an isosceles triangle pointing upward with the "center"
in the middle of the base, for the default parameter = `0.5`.

`Triangle(p)` constructs a triangle whose base is along the x-axis
going from (-p,0) to (1-p,0)` and with height=1.

The (width and height can be scaled when constructing an object).
"""
struct Triangle{T} <: AbstractShape2
    param::T # fraction in interval (0,1) 
    function Triangle{T}(param::T=0.5) where {T <: Real}
        0 < param < 1 || throw(ArgumentError("param=$param"))
        new{T}(param)
    end
end


#=
"""
    shape = Triangle(param=0.5)
"""
Triangle(param::T=0.5) where {T <: Real} = Triangle{T}(param)

function _trifun(x, y, param)
    return (x,y) -> (0 ≤ y ≤ 1) && (
end

function phantom(ob:Triangle)
    ob.shape == 1/2 || throw("todo")
    y
    return (x,y) -> (0 ≤ y ≤ 1) && (
end

# todo:
function phantom(ob::Object2d{Triangle})
    @show ob.shape.param
end



# tests
ob = Object(Triangle())


@show Triangle()
@show Triangle(3//4)
# @test_throws Triangle(2)
=#
