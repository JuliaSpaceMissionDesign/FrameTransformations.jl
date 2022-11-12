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

@testset "Function parse_iso" begin
    s = "2022-01-01T12:00"
    @test all(Tempo.parse_iso("2022-12-20") .== (2022, 12, 20, 0, 0, 0, 0))
    @test all(Tempo.parse_iso("2022-12-20T12") .== (2022, 12, 20, 12, 0, 0, 0))
    @test all(Tempo.parse_iso("2022-12-20T13:32") .== (2022, 12, 20, 13, 32, 0, 0))
    @test all(Tempo.parse_iso("2022-12-20T13:32:12") .== (2022, 12, 20, 13, 32, 12, 0))
    @test all(Tempo.parse_iso("2022-12-20T13:32:12.2132") .== (2022, 12, 20, 13, 32, 12, 0.2132))
    @test all(Tempo.parse_iso("2022-12-20T12:00:00.0 TDB") .== (2022, 12, 20, 12, 0, 0, 0))

    for _ in 1:5
        s, ry, rm, rd, rH, rM, rS, rF = generate_random_epoch()
        @test all(isapprox.(Tempo.parse_iso(s),(ry, rm, rd, rH, rM, rS, rF); rtol=1e-3)) 

        s, ry, rm, rd, rH, rM, rS, rF = generate_random_datetimeiso()
        @test all(isapprox.(Tempo.parse_iso(s),(ry, rm, rd, rH, rM, rS, rF); rtol=1e-3)) 
 
    end
end

@testset "Epoch" verbose=true begin
   
    @testset "String constructors" begin
        s, ry, rm, rd, rH, rM, rS, rF = generate_random_epoch()
        e = Epoch(s)
        dt = DateTime(e)
        @test all((ry, rm, rd, rH, rM, rS+rF) .≈ (year(dt), month(dt), 
            day(dt), hour(dt), minute(dt), second(dt)))
        @test Epoch("-0.5") ≈ Epoch("2000-01-01T00:00:00.0000 TDB")
        @test Epoch("0.5") ≈ Epoch("2000-01-02T00:00:00.0000 TDB")
        @test Epoch("JD 2400000.5") ≈ Epoch("1858-11-17T00:00:00.0000 TDB")
        @test Epoch("JD 2400000.5") ≈ Epoch("MJD 0.0")
    end
    
    @testset "DateTime and offset constructors" begin
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
    end

    # Based on Vallado "Fundamental of astrodynamics" page 196
    e = Epoch("2004-05-14T16:43:00 UTC")
    @test DateTime("2004-05-14T16:43:32.0000") ≈ DateTime(convert(TAI, e))
    @test DateTime("2004-05-14T16:44:04.1840") ≈ DateTime(convert(TT, e))
    @test DateTime("2004-05-14T16:44:04.1856") ≈ DateTime(convert(TDB, e))
    @test DateTime("2004-05-14T16:44:17.5255") ≈ DateTime(convert(TCB, e))
    @test DateTime("2004-05-14T16:44:04.7836") ≈ DateTime(convert(TCG, e))
end