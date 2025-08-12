module ImagePhantoms

const RealU = Number # Union{Real, Unitful.Length}

# core:
include("shape.jl")
include("object.jl")
include("jinc.jl")

# shapes:
include("ellipse.jl")
include("rect.jl")

end # module
