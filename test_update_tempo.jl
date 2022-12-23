const DJ2000 = 2451544.5

include("update-tempo/Tempo.jl")

using Test

### FIXME: move
### convert.jl 

@testset "Function isleapyear" begin
    for y in (2004, 2008, 2012, 2016)
        @test isleapyear(y)
    end

    for y in (2003, 2005, 2015, 2017)
        @test !isleapyear(y)
    end
end

@testset "Function find_dayinyear" begin 
    @test find_dayinyear(1, 24, false) == 24
    @test find_dayinyear(12, 30, false) == 364
end

@testset "Function hms2fd" begin 
    for _ in 1:10
        h = rand(0:23)
        m = rand(0:59)
        s = rand(0.0:59.0)
        @test hms2fd(h, m, s) ≈ ((s/60 + m)/60 + h)/24
    end

    @test_throws EpochConversionError hms2fd(24, 0, 0.0)
    @test_throws EpochConversionError hms2fd(23, 61, 0.0)
    @test_throws EpochConversionError hms2fd(20, 1, 61.0)
end

@testset "Function fd2hms" begin
    for _ in 1:10
        fd = rand()
        secinday = fd * 86400.0
        hours = Integer(secinday ÷ 3600)
        secinday -= 3600 * hours
        mins = Integer(secinday ÷ 60)
        secinday -= 60 * mins

        @test all(fd2hms(fd) .≈ (hours, mins, secinday))
        @test all(
            fd2hmsf(Float64(fd)) .≈ (hours, mins, secinday÷1, secinday-secinday÷1))
    end

    @test_throws EpochConversionError fd2hms(-0.1)
    @test_throws EpochConversionError fd2hms(1.5)
    @test_throws EpochConversionError fd2hmsf(-0.1)
    @test_throws EpochConversionError fd2hmsf(1.5)
end

@testset "Function cal2jd" begin
    Y, M, D = 2022, 6, 15
    h, m, s = 11, 51, 55.05
    
    @test sum(cal2jd(Y, M, D)) + hms2fd(h, m, s) ≈ 
        sum(calhms2jd(Y, M, D, h, m, s))

    @test_throws EpochConversionError cal2jd(rand(0:1580), rand(1:12), rand(1:28))
    @test_throws EpochConversionError cal2jd(rand(1600:2600), rand(1:12), 0)
    @test_throws EpochConversionError cal2jd(1850, 0, rand(1:28))
    @test_throws EpochConversionError cal2jd(2500, 13, 0)
    @test_throws EpochConversionError cal2jd(2150, 1, 32)
    @test_throws EpochConversionError cal2jd(1999, 6, 0)
    @test_throws MethodError cal2jd(1999, 13, 0.)
    @test_throws MethodError cal2jd(1999, 11., 0)
    
end

@testset "Function leapseconds" begin 
    for (y, m, l) in LEAP_TABLE[1:end]
        @test leapseconds(y, m) == l 
        @test leapseconds(y, m+1) == l 
        @test leapseconds(y, m-1) != l 
    end
end

@testset "Function utc2tai" begin
    for _ in 1:100
        Y, M, D = rand(1975:2015), rand(1:12), rand(1:28)
        h, m, s = rand(0:23), rand(0:59), rand(0.0:59.999)
        utc1, utc2 = calhms2jd(Y, M, D, h, m, s)
        tai1, tai2 = utc2tai(utc1, utc2)
        @test any(
            ((tai2 - utc2) * 86400 ≈ leapseconds(Y, M),
            (tai2 - utc2) * 86400 ≈ leapseconds(Y, M)+1))
    end
end

@testset "Function tai2utc" begin 

    for _ in 1:100
        Y, M, D = rand(1975:2015), rand(1:12), rand(1:28)
        h, m, s = rand(0:23), rand(0:59), rand(0.0:0.0001:59.999)
        utc1, utc2 = calhms2jd(Y, M, D, h, m, s)
        tai1, tai2 = utc2tai(utc1, utc2)
        u1, u2 = tai2utc(tai1, tai2)

        @test u1 ≈ utc1
        @test u2 ≈ utc2
    end

    # test limit case
    Y, M, D = 2008, 12, 31
    h, m, s = 23, 59, 59.99999999999999
    utc1, utc2 = calhms2jd(Y, M, D, h, m, s)
    tai1, tai2 = utc2tai(utc1, utc2)
    u1, u2 = tai2utc(tai1, tai2)

    @test u2 ≈ utc2

end

### FIXME: move
### parse.jl 

function _random_datetime()
    ry = rand(1800:2100)
    rm = rand(1:12)
    ly = isleapyear(ry)
    if ly && rm == 2
        rd = rand(1:29)
    elseif rm == 2
        rd = rand(1:28)
    else
        rd = rand(1:MTAB[rm])
    end
    rH = rand(0:23)
    rM = rand(0:59)
    rS = rand(0:59)
    rF = rand(0.:0.000000001:1.)

    return ry, rm, rd, rH, rM, rS, rF
end

function _random_datetime_isostr()
    ry, rm, rd, rH, rM, rS, rF = _random_datetime()
    rs = rS + rF
    return "$ry-$rm-$(rd)T$rH:$rM:$rs", ry, rm, rd, rH, rM, rS, rF
end


@testset "Function parse_iso" begin
    
    @test all(parse_iso("2021") .== (2021, 11, 1, 0, 0, 0, 0))
    @test all(parse_iso("2022-11") .== (2022, 11, 1, 0, 0, 0, 0))
    @test all(parse_iso("2022-12-20") .== (2022, 12, 20, 0, 0, 0, 0))
    @test all(parse_iso("2022-12-20T12") .== (2022, 12, 20, 12, 0, 0, 0))
    @test all(parse_iso("2022-12-20T13:32") .== (2022, 12, 20, 13, 32, 0, 0))
    @test all(parse_iso("2022-12-20T13:32:12") .== (2022, 12, 20, 13, 32, 12, 0))
    @test all(parse_iso("2022-12-20T13:32:12.2132") .== (2022, 12, 20, 13, 32, 12, 0.2132))
    @test all(parse_iso("2022-12-20T12:00:00.0 TDB") .== (2022, 12, 20, 12, 0, 0, 0))

    for _ in 1:100
        # s, ry, rm, rd, rH, rM, rS, rF = generate_random_epoch()
        # @test all(isapprox.(parse_iso(s),(ry, rm, rd, rH, rM, rS, rF); rtol=1e-3)) 

        s, ry, rm, rd, rH, rM, rS, rF = _random_datetime_isostr()
        @test all(isapprox.(parse_iso(s),(ry, rm, rd, rH, rM, rS, rF); rtol=1e-3)) 
 
    end
end


### FIXME: move
# offset.jl

@testset "UTC/TAI" begin 

    iytai = 2009; imtai = 1; idtai = 1;
    ihtai = 0; imintai = 0; sectai = 33.70
    tai1, tai2 = calhms2jd(iytai, imtai, idtai, ihtai, imintai, sectai)
    utc1, utc2 = tai2utc(tai1, tai2)
    @test sum(utc2tai(utc1, utc2)) == sum((tai1, tai2)) 

    utc1, utc2 = calhms2jd(2004, 5, 14, 10, 43, 0.0)
    tai1, tai2 = utc2tai(utc1, utc2)
    @test all(jd2calhms(tai1, tai2) .≈ (2004, 5, 14, 10, 43, 32.0))

    iytai = 2000; imtai = 1; idtai = 1;
    ihtai = 0; imintai = 0; sectai = rand(0.0:0.00001:59.999)
    tai1, tai2 = calhms2jd(iytai, imtai, idtai, ihtai, imintai, sectai)
    utc1, utc2 = tai2utc(tai1, tai2)
    @test offset_tai2utc(sectai) ≈ (utc2 - tai2) * 86400
    @test offset_utc2tai(utc2*86400) ≈ (tai2 - utc2) * 86400 

    utc1, utc2 = calhms2jd(2009, 1, 1, 0, 1, 0.7)
    tai1, tai2 = utc2tai(utc1, utc2)
    @test sum(tai2utc(tai1, tai2)) == sum((utc1, utc2))
    
end

@testset "TT/TAI" begin 
    @test offset_tai2tt(0.0) ≈ -offset_tt2tai(0.0)
    @test offset_tai2tt(0.0) == offset_tai2tt(rand(0.0:1e9))
    @test offset_tai2tt(0.0) == OFFSET_TAI_TT
end

@testset "TT/TCG" begin 
    @test offset_tt2tcg(0.0) ≈ -offset_tcg2tt(0.0)
end

@testset "TDB/TCB" begin 
    @test offset_tdb2tcb(0.0) + offset_tcb2tdb(0.0) < 1e-6
    @test offset_tdb2tcb(1e6) + offset_tcb2tdb(1e6) < 1e-6
end

@testset "TT/TDB" begin 
    @test offset_tt2tdb(0.0) ≈ -offset_tdb2tt(0.0)
    @test offset_tt2tdb(1e6) < 1/1000
    @test -offset_tdb2tt(1e6) < 1/1000
end

### FIXME: move
# scales.jl

function _one_offset(sec::Number)
    return sec+1.0 
end

function _mone_offset(sec::Number)
    return sec-1.0
end

S = TimeSystem{Float64}()

@timescale ATS 1 ATimeScale 
@timescale BTS 2 BTimeScale 
@timescale CTS 3 CTimeScale

add_timescale(S, ATS, _zero_offset)
add_timescale(S, BTS, _one_offset, parent=ATS, ftp=_mone_offset)
add_timescale(S, CTS, _one_offset, parent=BTS, ftp=_mone_offset)

@testset "TimeSystem" begin 

    @test apply_offsets(S, 0.0, ATS, CTS) ≈ 2.0 
    @test apply_offsets(S, 0.0, CTS, ATS) ≈ -2.0 
    @test apply_offsets(S, 0.0, ATS, BTS) ≈ 1.0 
    @test apply_offsets(S, 0.0, BTS, ATS) ≈ -1.0 

end
