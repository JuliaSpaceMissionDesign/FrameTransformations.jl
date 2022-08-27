@testset "Read generated file" begin
    gen = Utils.read_generated(joinpath(@__DIR__, "..", "assets", "test-generated.jl"))
    meta, body = parse_generated(gen)
    @test length(meta) == 2 
    @test length(body) == 2 
    @test meta[1] != meta[2]
end