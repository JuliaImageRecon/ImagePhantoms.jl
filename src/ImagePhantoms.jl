module ImagePhantoms

const RealU = Number # Union{Real, Unitful.Length}

# core:
include("shape.jl")
include("object.jl")
include("object3.jl")
include("jinc.jl")

# shapes:
include("ellipse.jl")
include("gauss2.jl")
include("rect.jl")
include("triangle.jl")

# 3d:
include("ellipsoid.jl")
include("gauss3.jl")
include("cuboid.jl")

# phantoms:
include("shepplogan.jl")
include("disk-phantom.jl")
include("focus-chart.jl")

# MRI SENSE:
include("mri-sense.jl")

end # module
