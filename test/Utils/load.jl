@testset "Load .tpc files" begin 
    data = load(TPC(joinpath(@__DIR__, "test.tpc")))
    @test data[399] == Dict(:gm => 3.9860043543609598E+05)
    @test length(keys(data[599])) == 7
    for (k, l, fv, lv) in zip(
        [:pole_ra, :pole_dec, :pm, :nut_prec_ra, :nut_prec_dec, :nut_prec_pm],
        [3, 3, 3, 15, 15, 15],
        [268.056595, 64.495303, 284.95, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.002150, 0.000926, 0.0],
    )
        @test length(data[599][k]) == l # test length
        @test data[599][k][1] == fv  # test first value
        @test data[599][k][end] == lv  # test last value
    end

    datas = load(
        [TPC(joinpath(@__DIR__, "test.tpc")), 
        TPC(joinpath(@__DIR__, "test2.tpc"))])

    @test datas[999][:gm] == -1.0
    @test datas[299][:radii][2] == 6051.8
    @test 5 in keys(datas)
end