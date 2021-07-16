# helper.jl
# utilities used in multiple test files

macro NOTinferred(ex) # flag where @inferred fails
    :($(esc(ex)))
end

# https://github.com/PainterQubits/Unitful.jl/issues/465
# reale = (x) -> (@assert x ≈ real(x); real(x))
myisreal = (z) -> maximum(abs ∘ imag, z) / maximum(abs, z) < sqrt(eps(eltype(real(z)))) * 10
reale = (z) -> (@assert myisreal(z); real(z))

# fft cannot handle units so this is a work-around
function myfft(x)
    u = unit(eltype(x))
    return fftshift(fft(fftshift(x) / u)) * u
end
