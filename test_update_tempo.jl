
const DJ2000 = 2451544.5

include("update-tempo/Tempo.jl")

using Test

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