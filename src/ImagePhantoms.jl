module ImagePhantoms

const RealU = Number # Union{Real, Unitful.Length}

# core:
include("shape.jl")
include("object.jl")
include("jinc.jl")

# shapes:
include("ellipse.jl")
include("gauss2.jl")
include("rect.jl")

# phantoms:
include("shepplogan.jl")

end # module
