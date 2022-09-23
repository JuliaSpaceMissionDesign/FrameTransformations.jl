@testset "convert.jl" begin 
    # Julian to calendar
    iy = 2008; im = 2; id = 29;
    ihour = 23; imin = 59; sec = 59.9;
    @test sum(Tempo.cal2jd(
        iy, im, id)) + Tempo.hms2fd(
            ihour, imin, sec)-0.5 ≈ 2454526.499999
    @test sum(Tempo.cal2jd(1858, 11, 17))-0.5 ≈ Tempo.DJM0

    fd = Tempo.hms2fd(ihour, imin, sec)
    jd1, jd2 = Tempo.cal2jd(iy, im, id)

    @test all(Tempo.jd2cal(jd1, jd2+fd-0.5) .≈ (iy, im, id, fd))
    @test all(Tempo.jd2calhms(
        Tempo.calhms2jd(
            iy, im, id, ihour, imin, sec)...) .≈ (iy, im, id, ihour, imin, sec))
    
    # output types 
    calhmsec = Tempo.jd2calhms(jd1, jd2+fd-0.5)
    @test all(typeof.(calhmsec) .== (Int, Int, Int, Int, Int, Float64))
    @test all(typeof.((jd1, jd1)) .== (Float64, Float64))
end

@testset "origin.jl" begin 
    iy = 1993; im = 2; id = 25;
    ihour = 12; imin = 53; sec = 15.0;

    # compute jd w.r.t J2000
    jd1, jd2 = Tempo.calhms2jd(iy, im, id, ihour, imin, sec)
    fd = Tempo.hms2fd(ihour, imin, sec)
    
    for origin in Tempo.EPOCH_ORIGIN_ACRONYMS
        # test parser
        @test Tempo.tryparse(Val(origin)) == eval(origin)

        jd1 -= Tempo.offset(Tempo.tryparse(Val(origin)))
        jd2 += Tempo.offset(Tempo.tryparse(Val(origin)))
        
        # transform back to date  
        @test all(Tempo.jd2cal(jd1, jd2) .≈ (iy, im, id, fd))
    end
end




@testset "scales.jl" begin
    # check if all declared scales are inserted in the graph
    scales = collect(keys(TIMESCALES.nodes))
    for scale in Tempo.TIMESCALES_ACRONYMS
        @test tryparse(Val(scale)) in scales 
    end

    iytai = 2009; imtai = 1; idtai = 1;
    ihtai = 0; imintai = 0; sectai = 33.70
    tai1, tai2 = Tempo.calhms2jd(iytai, imtai, idtai, ihtai, imintai, sectai)
    utc1, utc2 = Tempo.tai2utc(tai1, tai2)
    @test sum(Tempo.utc2tai(utc1, utc2)) == sum((tai1, tai2)) 

    utc1, utc2 = Tempo.calhms2jd(2009, 1, 1, 0, 1, 0.7)
    tai1, tai2 = Tempo.utc2tai(utc1, utc2)
    @test sum(Tempo.tai2utc(tai1, tai2)) == sum((utc1, utc2))

    utc1, utc2 = Tempo.calhms2jd(2004, 5, 14, 10, 43, 0.0)
    tai1, tai2 = Tempo.utc2tai(utc1, utc2)
    @test all(Tempo.jd2calhms(tai1, tai2) .≈ (2004, 5, 14, 10, 43, 32.0))

    @testset "offset.jl" begin 

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
end

@testset "parse.jl" begin 
    s = "2022-01-01T12:00"
    @test all(Tempo.parse_iso("2022-12-20") .== (2022, 12, 20, 0, 0, 0, 0))
    @test all(Tempo.parse_iso("2022-12-20T12") .== (2022, 12, 20, 12, 0, 0, 0))
    @test all(Tempo.parse_iso("2022-12-20T13:32") .== (2022, 12, 20, 13, 32, 0, 0))
    @test all(Tempo.parse_iso("2022-12-20T13:32:12") .== (2022, 12, 20, 13, 32, 12, 0))
    @test all(Tempo.parse_iso("2022-12-20T13:32:12.2132") .== (2022, 12, 20, 13, 32, 12, 0.2132))
    @test all(Tempo.parse_iso("2022-12-20T12:00:00.0 TDB") .== (2022, 12, 20, 12, 0, 0, 0))
end

@testset "datetime.jl" begin 
    @testset "Date" begin
        d = Date(2000, 1, 1)

        @test year(d) == 2000
        @test month(d) == 1 
        @test day(d) == 1
        @test isleapyear(d) 
        @test find_dayinyear(d) == 1 

        @test d == Date(0)
        @test d == Date(2000, 1)
        @test j2000(d) == 0

        @test Tempo.find_month(32, true) == 2
        @test Tempo.find_day(32, 2, true) == 1
        @test d + 1 == Date(1)
        @test Date(1) - 1 == d 
    end
    @testset "Time" begin
        t0 = Time(12, 0, 0.0)
        @test t0 == Time(86400÷2, 0.0)

        t = Time(12, 1, 15.300300300)
        @test Tempo.millisecond(t) == 300
        @test Tempo.microsecond(t) == 300
        @test Tempo.nanosecond(t) == 300
        @test Tempo.second(Float64, t) == 15.300300300
        @test Tempo.second(Int64, t) == 15
        @test Tempo.fraction_of_second(t) ≈ 0.300300300
    end
end


# function Epoch(ep::Epoch{S1}, ::S2) where {S1<:TimeScale, S2<:TimeScale}
#     second, fraction, error = apply_offset(ep.second, ep.fraction, ep.error, S1(), S2())
#     Epoch{S2}(second, fraction, error)
# end

# Epoch(ep::Epoch{S}, ::S) where {S<:TimeScale} = ep

# j2000seconds(e::Epoch) = value(Epoch(e, TDB))
# j2000(e::Epoch) = j2000seconds(e)/86400.0
# j2000centuries(e::Epoch) = j2000seconds(e)/(86400.0*365.25*100.0)
