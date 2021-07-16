macro NOTinferred(ex) # flag where @inferred fails
    :($(esc(ex)))
end
