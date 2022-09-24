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

function generate_random_datetime()
    ry = rand(1800:2100)
    rm = rand(1:12)
    ly = Tempo.isleapyear(ry)
    if ly && rm == 2
        rd = rand(1:29)
    elseif rm == 2
        rd = rand(1:28)
    else
        rd = rand(1:Tempo.MTAB[rm])
    end
    rH = rand(0:23)
    rM = rand(0:59)
    rS = rand(0:59)
    rF = rand(0.:0.000000001:1.)

    return ry, rm, rd, rH, rM, rS, rF
end

function generate_random_datetimeiso()
    ry, rm, rd, rH, rM, rS, rF = generate_random_datetime()
    rs = rS + rF
    return "$ry-$rm-$(rd)T$rH:$rM:$rs", ry, rm, rd, rH, rM, rS, rF
end

function generate_random_epoch()
    ry, rm, rd, rH, rM, rS, rF = generate_random_datetime()
    rs = rS + rF
    nscales = length(Tempo.TIMESCALES_ACRONYMS)
    rscale = rand(1:nscales)
    return "$ry-$rm-$(rd)T$rH:$rM:$rs $(String(Tempo.TIMESCALES_ACRONYMS[rscale]))", ry, rm, rd, rH, rM, rS, rF
end

@testset "parse.jl" begin 
    s = "2022-01-01T12:00"
    @test all(Tempo.parse_iso("2022-12-20") .== (2022, 12, 20, 0, 0, 0, 0))
    @test all(Tempo.parse_iso("2022-12-20T12") .== (2022, 12, 20, 12, 0, 0, 0))
    @test all(Tempo.parse_iso("2022-12-20T13:32") .== (2022, 12, 20, 13, 32, 0, 0))
    @test all(Tempo.parse_iso("2022-12-20T13:32:12") .== (2022, 12, 20, 13, 32, 12, 0))
    @test all(Tempo.parse_iso("2022-12-20T13:32:12.2132") .== (2022, 12, 20, 13, 32, 12, 0.2132))
    @test all(Tempo.parse_iso("2022-12-20T12:00:00.0 TDB") .== (2022, 12, 20, 12, 0, 0, 0))

    s, ry, rm, rd, rH, rM, rS, rF = generate_random_epoch()
    @test all(isapprox.(Tempo.parse_iso(s),(ry, rm, rd, rH, rM, rS, rF); rtol=1e-3)) 

    s, ry, rm, rd, rH, rM, rS, rF = generate_random_datetimeiso()
    @test all(isapprox.(Tempo.parse_iso(s),(ry, rm, rd, rH, rM, rS, rF); rtol=1e-3)) 
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

    @testset "DateTime" begin
        d = Date(2000, 1, 1)
        t0 = Time(12, 0, 0.0)
        J2000 = DateTime(d, 0.0)
        @test t0 == Time(J2000)
        @test d == Date(J2000)

        s = "2000-01-01T12:00:00.0"
        @test J2000 == DateTime(s)
        @test J2000 == DateTime(0.0)

        D2 = DateTime(Tempo.DAY2SEC/3)
        @test Tempo.j2000(D2) ≈ 1/3
        @test Tempo.j2000s(D2) ≈ 1/3 * Tempo.DAY2SEC
        @test Tempo.j2000c(D2) ≈ 1/3 / Tempo.CENTURY2DAY

        ry, rm, rd, rH, rM, rS, rF = generate_random_datetime()
        dt = DateTime(ry, rm, rd, rH, rM, rS, rF)
        @test year(dt) == ry 
        @test month(dt) == rm
        @test day(dt) == rd
        @test hour(dt) == rH
        @test minute(dt) == rM
        @test second(dt) == rS + rF

        dtr = rand(0.0:1.0:365.25*86400.0)
        dt2 = dt + dtr
        @test Tempo.j2000s(dt2) - Tempo.j2000s(dt) ≈ dtr
    end
    @testset "Epoch" begin
        s, ry, rm, rd, rH, rM, rS, rF = generate_random_epoch()
        e = Epoch(s)
        dt = DateTime(e)

        @test all((ry, rm, rd, rH, rM, rS+rF) .≈ (year(dt), month(dt), 
            day(dt), hour(dt), minute(dt), second(dt)))
        @test Epoch("-0.5") ≈ Epoch("2000-01-01T00:00:00.0000 TDB")
        @test Epoch("0.5") ≈ Epoch("2000-01-02T00:00:00.0000 TDB")
        @test Epoch("JD 2400000.5") ≈ Epoch("1858-11-17T00:00:00.0000 TDB")
        @test Epoch("JD 2400000.5") ≈ Epoch("MJD 0.0")
        
        s, ry, rm, rd, rH, rM, rS, rF = generate_random_epoch()
        e = Epoch(s)
        dt = DateTime(e)
        @test DateTime(e) ≈ dt

        rn = rand(0:10000)
        @test value(e + rn) ≈ value(e) + rn
        rn = rand(0:10000)
        @test value(e - rn) ≈ value(e) - rn

        rn0 = rand(-2000:2000)
        rn1 = rand(-1000:1000)
        @test Epoch("$rn0") - Epoch("$rn1") ≈ (rn0 - rn1)*Tempo.DAY2SEC
        @test all(
            collect(Epoch("0"):86400.0:Epoch("2")) .== [Epoch("0"), Epoch("1"), Epoch("2")])

        # Based on Vallado "Fundamental of astrodynamics" page 196
        e = Epoch("2004-05-14T16:43:00 UTC")
        @test DateTime("2004-05-14T16:43:32.0000") ≈ DateTime(convert(TAI, e))
        @test DateTime("2004-05-14T16:44:04.1840") ≈ DateTime(convert(TT, e))
        @test DateTime("2004-05-14T16:44:04.1856") ≈ DateTime(convert(TDB, e))
        @test DateTime("2004-05-14T16:44:17.5255") ≈ DateTime(convert(TCB, e))
        @test DateTime("2004-05-14T16:44:04.7836") ≈ DateTime(convert(TCG, e))
    end
end