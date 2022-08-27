@testset "Generate unique fileids" begin
    fp1 = joinpath(@__DIR__, "..", "assets", "test.tpc")
    fp2 = joinpath(@__DIR__, "..", "assets", "test2.tpc")

    fid1 = Basic.Utils.fileid(fp1)
    fid2 = Basic.Utils.fileid(fp2)
    fid1_ = Basic.Utils.fileid(CONFIG(fp1))
    @test fid1 == fid1_ 
    @test fid2 != fid1
end