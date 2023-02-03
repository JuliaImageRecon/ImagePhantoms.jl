module ImagePhantoms

const RealU = Number # Union{Real, Unitful.Length}

# core:
include("shape.jl")
include("rotate3.jl")
include("object.jl")
include("object2.jl")
include("object3.jl")
include("jinc.jl")

# 2d:
include("dirac2.jl")
include("ellipse.jl")
include("gauss2.jl")
include("rect.jl")
include("triangle.jl")

# 3d:
include("cone.jl")
include("cuboid.jl")
include("cylinder.jl")
include("dirac3.jl")
include("ellipsoid.jl")
include("gauss3.jl")

# phantoms:
include("shepplogan.jl")
include("disk-phantom.jl")
include("focus-chart.jl")

# MRI SENSE:
include("mri-sense.jl")

end # module
