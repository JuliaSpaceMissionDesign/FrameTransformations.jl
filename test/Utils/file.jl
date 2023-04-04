@testset "Generate unique fileids" begin
    fp1 = joinpath(@__DIR__, "..", "assets", "gms.tpc")
    fp2 = joinpath(@__DIR__, "..", "assets", "orient.tpc")
    fid1 = Utils.fileid(fp1)
    fid2 = Utils.fileid(fp2)
    @test fid2 != fid1
end

@testset "Function filepath" begin
    file = TPC(joinpath(@__DIR__, "..", "assets", "gms.tpc"))
    @test Utils.filepath(file) == joinpath(@__DIR__, "..", "assets", "gms.tpc")
end
