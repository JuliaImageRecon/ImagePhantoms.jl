# helper.jl
# utilities used in multiple test files

#import FFTW: fft

#=
macro NOTinferred(ex) # flag where @inferred fails
    :($(esc(ex)))
end
=#


if !@isdefined myfft
    # fft cannot handle units so this is a work-around
    function myfft(x::AbstractArray{<:Number})
        u = oneunit(eltype(x))
        return fftshift(fft(ifftshift(x) / u)) * u
    end
end
