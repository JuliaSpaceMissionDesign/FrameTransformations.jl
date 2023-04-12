@testset "Offset functions" verbose = true begin
    @testset "UTC/TAI" begin
        iytai = 2009
        imtai = 1
        idtai = 1
        ihtai = 0
        imintai = 0
        sectai = 34.0
        tai1, tai2 = Tempo.calhms2jd(iytai, imtai, idtai, ihtai, imintai, sectai)
        utc1, utc2 = Tempo.tai2utc(tai1, tai2)
        @test sum(Tempo.utc2tai(utc1, utc2)) == sum((tai1, tai2))

        utc1, utc2 = Tempo.calhms2jd(2004, 5, 14, 10, 43, 0.0)
        tai1, tai2 = Tempo.utc2tai(utc1, utc2)
        @test all(Tempo.jd2calhms(tai1, tai2) .≈ (2004, 5, 14, 10, 43, 32.0))

        iytai = 2000
        imtai = 1
        idtai = 1
        ihtai = 0
        imintai = 0
        sectai = rand(0.0:0.00001:59.999)
        tai1, tai2 = Tempo.calhms2jd(iytai, imtai, idtai, ihtai, imintai, sectai)
        utc1, utc2 = Tempo.tai2utc(tai1, tai2)
        @test Tempo.offset_tai2utc(sectai) ≈ (utc2 - tai2) * 86400
        @test Tempo.offset_utc2tai(utc2 * 86400) ≈ (tai2 - utc2) * 86400

        utc1, utc2 = Tempo.calhms2jd(2009, 1, 1, 0, 1, 0.7)
        tai1, tai2 = Tempo.utc2tai(utc1, utc2)
        @test sum(Tempo.tai2utc(tai1, tai2)) == sum((utc1, utc2))
    end

    @testset "TT/TAI" begin
        @test Tempo.offset_tai2tt(0.0) ≈ -Tempo.offset_tt2tai(0.0)
        @test Tempo.offset_tai2tt(0.0) == Tempo.offset_tai2tt(rand(0.0:1e9))
        @test Tempo.offset_tai2tt(0.0) == Tempo.OFFSET_TAI_TT
    end

    @testset "TT/TCG" begin
        @test Tempo.offset_tt2tcg(0.0) ≈ -Tempo.offset_tcg2tt(0.0)
    end

    @testset "TDB/TCB" begin
        @test Tempo.offset_tdb2tcb(0.0) + Tempo.offset_tcb2tdb(0.0) < 1e-6
        @test Tempo.offset_tdb2tcb(1e6) + Tempo.offset_tcb2tdb(1e6) < 1e-6
    end

    @testset "TT/TDB" begin
        @test Tempo.offset_tt2tdb(0.0) ≈ -Tempo.offset_tdb2tt(0.0)
        @test Tempo.offset_tt2tdb(1e6) < 1 / 1000
        @test -Tempo.offset_tdb2tt(1e6) < 1 / 1000
    end
end
