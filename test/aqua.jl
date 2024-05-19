using ImagePhantoms: ImagePhantoms
import Aqua
using Test: @testset

@testset "aqua" begin
    Aqua.test_all(ImagePhantoms)
end
