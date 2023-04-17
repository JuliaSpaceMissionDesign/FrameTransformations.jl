
@testset "Function parse_iso" begin
    @test all(parse_iso("2021") .== (2021, 1, 1, 0, 0, 0, 0))
    @test all(parse_iso("2022-11") .== (2022, 11, 1, 0, 0, 0, 0))
    @test all(parse_iso("2022-12-20") .== (2022, 12, 20, 0, 0, 0, 0))
    @test all(parse_iso("2022-12-20T12") .== (2022, 12, 20, 12, 0, 0, 0))
    @test all(parse_iso("2022-12-20T13:32") .== (2022, 12, 20, 13, 32, 0, 0))
    @test all(parse_iso("2022-12-20T13:32:12") .== (2022, 12, 20, 13, 32, 12, 0))
    @test all(parse_iso("2022-12-20T13:32:12.2132") .== (2022, 12, 20, 13, 32, 12, 0.2132))
    @test all(parse_iso("2022-12-20T12:00:00.0 TDB") .== (2022, 12, 20, 12, 0, 0, 0))

    for _ in 1:10
        # s, ry, rm, rd, rH, rM, rS, rF = generate_random_epoch()
        # @test all(isapprox.(parse_iso(s),(ry, rm, rd, rH, rM, rS, rF); rtol=1e-3)) 

        s, ry, rm, rd, rH, rM, rS, rF = _random_datetime_isostr()
        @test all(isapprox.(parse_iso(s), (ry, rm, rd, rH, rM, rS, rF); rtol=1e-3))
    end
end
