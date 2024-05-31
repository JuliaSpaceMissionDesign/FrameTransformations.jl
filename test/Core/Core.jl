using Test 
using SafeTestsets
 
@testset "Core" verbose=true begin
    @safetestset "Rotations" begin include("rotation.jl") end
    @safetestset "Transform: points" begin include("transform_points.jl") end
    @safetestset "Transform: axes" begin include("transform_axes.jl") end
    @safetestset "Graph: points" begin include("points.jl") end
    @safetestset "Graph: frames" begin include("graph.jl") end
    @safetestset "Graph: directions" begin include("directions.jl") end
end;