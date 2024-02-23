using Test 

@testset "FrameTransformations" verbose=true begin
    for group in [:Core, :Definitions]
        include("$group/$group.jl")
    end
end;
