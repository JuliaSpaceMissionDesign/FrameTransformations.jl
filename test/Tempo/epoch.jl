
@testset "Epoch" verbose = true begin
    @testset "String constructors" begin
        s, ry, rm, rd, rH, rM, rS, rF = _random_epoch()
        e = Epoch(s)
        dt = DateTime(e)
        @test Epoch("-0.5") ≈ Epoch("2000-01-01T00:00:00.0000 TDB")
        @test Epoch("0.5") ≈ Epoch("2000-01-02T00:00:00.0000 TDB")
        @test Epoch("JD 2400000.5") ≈ Epoch("1858-11-17T00:00:00.0000 TDB")
        @test Epoch("JD 2400000.5") ≈ Epoch("MJD 0.0")
    end

    @testset "DateTime and offset constructors" begin
        s, ry, rm, rd, rH, rM, rS, rF = _random_epoch()
        e = Epoch(s)
        dt = DateTime(e)
        @test DateTime(e) ≈ dt

        rn = rand(0:10000)
        @test value(e + rn) ≈ value(e) + rn
        rn = rand(0:10000)
        @test value(e - rn) ≈ value(e) - rn

        rn0 = rand(-2000:2000)
        rn1 = rand(-1000:1000)
        @test Epoch("$rn0") - Epoch("$rn1") ≈ (rn0 - rn1) * Tempo.DAY2SEC
        @test all(
            collect(Epoch("0"):86400.0:Epoch("2")) .== [Epoch("0"), Epoch("1"), Epoch("2")]
        )

        @test Epoch(0.0, TDB) ≈ Epoch("2000-01-01T12:00:00.0000 TDB")
    end

    # Based on Vallado "Fundamental of astrodynamics" page 196
    e = Epoch("2004-05-14T16:43:00 UTC")
    @test DateTime("2004-05-14T16:43:32.0000") ≈ DateTime(convert(TAI, e))
    @test DateTime("2004-05-14T16:44:04.1840") ≈ DateTime(convert(TT, e))
    @test DateTime("2004-05-14T16:44:04.1856") ≈ DateTime(convert(TDB, e))
    @test DateTime("2004-05-14T16:44:17.5255") ≈ DateTime(convert(TCB, e))
    @test DateTime("2004-05-14T16:44:04.7836") ≈ DateTime(convert(TCG, e))
end
