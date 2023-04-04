function _one_offset(sec::Number)
    return 1.0
end

function _mone_offset(sec::Number)
    return -1.0
end

import Basic.Tempo: timescale_id, timescale_name

S = Tempo.TimeSystem{Float64}()
@timescale ATS 1 ATimeScale
@timescale BTS 2 BTimeScale
@timescale CTS 3 CTimeScale

Tempo.add_timescale(S, ATS, Tempo._zero_offset)
Tempo.add_timescale(S, BTS, _one_offset; parent=ATS, ftp=_mone_offset)
Tempo.add_timescale(S, CTS, _one_offset; parent=BTS, ftp=_mone_offset)

@testset "TimeSystem" verbose = true begin
    @test Tempo.apply_offsets(S, 0.0, ATS, CTS) ≈ 2.0
    @test Tempo.apply_offsets(S, 0.0, CTS, ATS) ≈ -2.0
    @test Tempo.apply_offsets(S, 0.0, ATS, BTS) ≈ 1.0
    @test Tempo.apply_offsets(S, 0.0, BTS, ATS) ≈ -1.0

    @testset "apply_offsets" begin
        D2S = 86400.0

        utc1, utc2 = Tempo.calhms2jd(2022, 1, 1, 12, 0, 0.0)
        tai1, tai2 = Tempo.utc2tai(utc1, utc2)
        @test (tai2 - utc2) * D2S ≈ Tempo.leapseconds(utc2)

        # Seconds since j2000
        j2000s_utc = utc2 * D2S

        # TAI
        @test Tempo.apply_offsets(TIMESCALES, j2000s_utc, UTC, TAI) - j2000s_utc ≈
            Tempo.leapseconds(utc2)
        @test Tempo.apply_offsets(TIMESCALES, j2000s_utc, UTC, TAI) - j2000s_utc ≈ 37.0

        # Same scale (no offset applied)
        for s in Tempo.TIMESCALES_ACRONYMS
            @eval begin
                @test Tempo.apply_offsets(TIMESCALES, 0.0, $s, $s) == 0.0
            end
        end

        # TT 
        @test Tempo.apply_offsets(TIMESCALES, 0.0, TAI, TT) == 32.184
        @test Tempo.apply_offsets(TIMESCALES, j2000s_utc, UTC, TT) - j2000s_utc ≈ 69.184
        @test Tempo.apply_offsets(
            TIMESCALES, Tempo.apply_offsets(TIMESCALES, j2000s_utc, UTC, TAI), TAI, TT
        ) ≈ Tempo.apply_offsets(TIMESCALES, j2000s_utc, UTC, TT)

        # Example from Vallado 
        jd2000, dutc = Tempo.calhms2jd(2004, 5, 14, 10, 43, 0.0)
        sutc = dutc * D2S

        # UTC - TAI
        stai = Tempo.apply_offsets(TIMESCALES, sutc, UTC, TAI)
        dtai = stai / D2S
        @test all(Tempo.jd2calhms(jd2000, dtai) .≈ (2004, 5, 14, 10, 43, 32.0))

        # UTC - TT 
        stt = Tempo.apply_offsets(TIMESCALES, sutc, UTC, TT)
        dtt = stt / D2S
        @test all(Tempo.jd2calhms(jd2000, dtt) .≈ (2004, 5, 14, 10, 44, 4.1840))
        @test Tempo.jd2calhms(
            jd2000, Tempo.apply_offsets(TIMESCALES, stai, TAI, TT) / D2S
        ) == Tempo.jd2calhms(jd2000, dtt)

        # UTC - TCB 
        stcb = Tempo.apply_offsets(TIMESCALES, sutc, UTC, TCB)
        dtcb = stcb / D2S
        @test all(
            isapprox.(
                Tempo.jd2calhms(jd2000, dtcb), (2004, 5, 14, 10, 44, 17.5255); rtol=1e-2
            ),
        )

        # UTC - TDB 
        stdb = Tempo.apply_offsets(TIMESCALES, sutc, UTC, TDB)
        dtdb = stdb / D2S
        @test all(
            isapprox.(
                Tempo.jd2calhms(jd2000, dtdb), (2004, 5, 14, 10, 44, 4.1856); rtol=1e-2
            ),
        )
    end
end
