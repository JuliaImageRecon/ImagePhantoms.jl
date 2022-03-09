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
include("triangle.jl")

include("cuboid.jl")

# phantoms:
include("shepplogan.jl")
include("disk-phantom.jl")
include("focus-chart.jl")

# MRI SENSE:
include("mri-sense.jl")

end # module
