#=
disk-phantom.jl
generate random disk phantoms
2019-07-03, Jeff Fessler, University of Michigan
=#

export disk_phantom_params

using Random: seed!


"""
    params = disk_phantom_params( ; ...)

Generate `ndisk` ellipse phantom parameter 6-tuples
for a head-sized disk plus many disks within it,
designed so that the disks have some minimum separation `minsep`
to avoid overlap and to simplify patch-based model fitting.

# Options
- `fov::Real = 240` image field of view in mm
- `rhead::Function = () -> 100` background radius for "head" [mm]
- `muhead::Function = () -> 1000` "μ" (intensity) value for background head disk
- `rdisk::Function = () -> 10 + 10 * rand()` random disk radii [10,20]
- `mudisk::Function = () -> 100 + 200 * rand()` "μ" values for disks [100,300]
- `ndisk::Function = () -> 10` # of random disks
- `minsep::Real = 8` minimum disk separation in mm
- `maxtry::Int = 500` give up on adding more disks if this is reached
- `warn::Bool = false` warn if maxtry reached?
- `seed::Int = 0` if nonzero then use this seed

The function options can be replaced with rand() for other Distributions.
"""
function disk_phantom_params( ; kwargs...)
    params = _disk_phantom_params( ; kwargs...)

    out = ellipse_parameters_uscale(params, (1,1), 1, 1, 1)
    return out
end


function _disk_phantom_params( ;
    fov::RealU = 240,
#   rhead::RealU = 100, # background radius for "head" [mm]
    rhead::Function = () -> 100, # background radius for "head" [mm]
    muhead::Function = () -> 1000, # "μ" value for background head
#   rmin::RealU = 10, # min radius for disks
#   rmax::RealU = 20, # max radius for disks
    rdisk::Function = () -> 10 + 10 * rand(), # disk radius [10,20]
#   mumin::RealU = 100, # range of "μ" values for disks
#   mumax::RealU = 300,
    mudisk::Function = () -> 100 + 200 * rand(), # "μ" values for disks [100,300]
    ndisk::Function = () -> rand(10:10), # of disks
    minsep::RealU = 8, # minimum separation in mm
    # radial position of disk centers, uses rhead and rmin:
    cdisk::Function = () -> sqrt(rand()) * (100 - 10 - minsep),
    maxtry::Int = 500, # give up on adding more disks if this is reached
    warn::Bool = false, # warn if maxtry reached?
    seed::Int = 0,
)

    ndisk = ndisk()
    rhead = rhead()
    μhead = muhead()

    params = zeros(Float32, ndisk+1, 6)
    params[end,:] = [0, 0, rhead, rhead, 0, μhead]

#   randu = (a, b ; f = x->x) -> a + f(rand()) * (b - a)
    seper = (a,b,r1,x,y,r2) -> sqrt(abs2(x - a) + abs2(y - b)) - r1 - r2

    (seed != 0) && seed!(seed)

    for id in 1:ndisk
        trial = 0
        for ii in 1:maxtry
#           rc = randu(0, rhead - rmin - minsep, f = x->x^0.5)
            rc = cdisk()
            phi = rand() * 2π
            rad = rdisk()
            val = mudisk()
            trial = [rc*cos(phi), rc*sin(phi), rad, rad, 0, val]

            id == 1 && break # first one is always fine

            # see if trial is too close to another one
            c = 1:(id-1) # check these
            sep = seper.(params[c,1], params[c,2], params[c,3],
                trial[1], trial[2], trial[3])

            (minimum(sep) > minsep) && (rhead - (rc + rad) > minsep) && break
            ii == maxtry && warn && @warn("need more tries")
            ii == maxtry && return params[[1:id; end],:]
        end

#=
        if id > 1
            c = 1:(id-1) # check these
            sep = seper(params[c,1], params[c,2], params[c,3],
                trial[1], trial[2], trial[3])
            @show minimum(sep)
        end
        @show id, trial
=#

        params[id,:] = trial
    end

    return params
end
