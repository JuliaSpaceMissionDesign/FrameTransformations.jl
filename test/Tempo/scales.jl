@testset "Check scales connection" begin
    # check if all declared scales are inserted in the scales graph
    scales = collect(keys(TIMESCALES.nodes))
    for scale in Tempo.TIMESCALES_ACRONYMS
        @test tryparse(Val(scale)) in scales 
    end
end

@testset "UTC to TAI" begin 
    iytai = 2009; imtai = 1; idtai = 1;
    ihtai = 0; imintai = 0; sectai = 33.70
    tai1, tai2 = Tempo.calhms2jd(iytai, imtai, idtai, ihtai, imintai, sectai)
    utc1, utc2 = Tempo.tai2utc(tai1, tai2)
    @test sum(Tempo.utc2tai(utc1, utc2)) == sum((tai1, tai2)) 

    utc1, utc2 = Tempo.calhms2jd(2004, 5, 14, 10, 43, 0.0)
    tai1, tai2 = Tempo.utc2tai(utc1, utc2)
    @test all(Tempo.jd2calhms(tai1, tai2) .≈ (2004, 5, 14, 10, 43, 32.0))
end

@testset "TAI to UTC" begin 
    utc1, utc2 = Tempo.calhms2jd(2009, 1, 1, 0, 1, 0.7)
    tai1, tai2 = Tempo.utc2tai(utc1, utc2)
    @test sum(Tempo.tai2utc(tai1, tai2)) == sum((utc1, utc2))
end

@testset "Scales offsets" begin 
    function apply_offsets(second, from, to)
        sint = floor(Int64, second)
        sfrc = second - sint 
        sum(Tempo.apply_offsets(sint, sfrc, from, to))
    end 

    D2S = 86400.0

    # test conversions 
    # Based on Vallado "Fundamental of astrodynamics" page 196
    jd2000, dutc = Tempo.calhms2jd(2004, 5, 14, 10, 43, 0.0)
    sutc = dutc * D2S

    stai = apply_offsets(sutc, UTC, TAI)
    dtai = stai / D2S
    @test all(Tempo.jd2calhms(jd2000, dtai) .≈ (2004, 5, 14, 10, 43, 32.0))
    
    stt = apply_offsets(sutc, UTC, TT)
    dtt = stt / D2S
    @test all(Tempo.jd2calhms(jd2000, dtt) .≈ (2004, 5, 14, 10, 44, 4.1840))
    @test Tempo.jd2calhms(jd2000, 
        apply_offsets(stai, TAI, TT)/D2S) == Tempo.jd2calhms(jd2000, dtt)

    stcb = apply_offsets(sutc, UTC, TCB)
    dtcb = stcb / D2S
    @test all(isapprox.(Tempo.jd2calhms(jd2000, dtcb),
        (2004, 5, 14, 10, 44, 17.5255); rtol=1e-1))

    stdb = apply_offsets(sutc, UTC, TDB)
    dtdb = stdb / D2S
    @test all(isapprox.(Tempo.jd2calhms(jd2000, dtdb),
        (2004, 5, 14, 10, 44, 4.1856); rtol=1e-1))
end