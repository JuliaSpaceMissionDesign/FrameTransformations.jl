@testset "Function isleapyear" begin
    for y in (2004, 2008, 2012, 2016)
        @test Tempo.isleapyear(y)
    end

    for y in (2003, 2005, 2015, 2017)
        @test !Tempo.isleapyear(y)
    end
end

@testset "Function find_dayinyear" begin 
    @test Tempo.find_dayinyear(1, 24, false) == 24
    @test Tempo.find_dayinyear(12, 30, false) == 364
end

@testset "Function hms2fd" begin 
    for _ in 1:10
        h = rand(0:23)
        m = rand(0:59)
        s = rand(0.0:59.0)
        @test Tempo.hms2fd(h, m, s) ≈ ((s/60 + m)/60 + h)/24
    end

    @test_throws ErrorException Tempo.hms2fd(24, 0, 0.0)
    @test_throws ErrorException Tempo.hms2fd(23, 61, 0.0)
    @test_throws ErrorException Tempo.hms2fd(20, 1, 61.0)
end

@testset "Function fd2hms" begin
    for _ in 1:10
        fd = rand()
        secinday = fd * 86400.0
        hours = Integer(secinday ÷ 3600)
        secinday -= 3600 * hours
        mins = Integer(secinday ÷ 60)
        secinday -= 60 * mins

        @test all(Tempo.fd2hms(fd) .≈ (hours, mins, secinday))
        @test all(
            Tempo.fd2hmsf(Float64(fd)) .≈ (hours, mins, secinday÷1, secinday-secinday÷1))
    end

    @test_throws ArgumentError Tempo.fd2hms(-0.1)
    @test_throws ArgumentError Tempo.fd2hms(1.5)
    @test_throws ArgumentError Tempo.fd2hmsf(-0.1)
    @test_throws ArgumentError Tempo.fd2hmsf(1.5)
end

@testset "Function cal2jd" begin
    iy = 2008; im = 2; id = 29;
    ihour = 23; imin = 59; sec = 59.9;
    @test sum(Tempo.cal2jd(
        iy, im, id)) + Tempo.hms2fd(
            ihour, imin, sec)-0.5 ≈ 2454526.499999

    @test sum(Tempo.cal2jd(1858, 11, 17))-0.5 ≈ Tempo.DJM0

    @test_throws ErrorException Tempo.cal2jd(1000, 0, 0)
    @test_throws ErrorException Tempo.cal2jd(1999, 0, 0)
    @test_throws ErrorException Tempo.cal2jd(1999, 13, 0)
    @test_throws MethodError Tempo.cal2jd(1999, 13, 0.)
    @test_throws MethodError Tempo.cal2jd(1999, 11., 0)
    @test_throws ErrorException Tempo.cal2jd(1999, 1, 32)
    @test_throws ErrorException Tempo.cal2jd(1999, 2, 30)
    @test_throws ErrorException Tempo.cal2jd(1999, 6, 0)
end

@testset "Funcion jd2cal" begin 
    iy = 2008; im = 2; id = 29;
    ihour = 23; imin = 59; sec = 59.9;
    jd1, jd2 = Tempo.cal2jd(iy, im, id)
    fd = Tempo.hms2fd(ihour, imin, sec)
    @test all(Tempo.jd2cal(jd1, jd2+fd-0.5) .≈ (iy, im, id, fd))
    @test all(Tempo.jd2calhms(
        Tempo.calhms2jd(
            iy, im, id, ihour, imin, sec)...) .≈ (iy, im, id, ihour, imin, sec))

    @test_throws ErrorException Tempo.jd2cal(0, 1e12)
    @test_throws ErrorException Tempo.jd2cal(0, -70000.0)

    # default return fraction of a day of 0.5 in case integer epochs given
    @test all(Tempo.jd2cal(Tempo.DJ2000, 366.0) .≈ (2001, 1, 1, 0.5))
    Y, M, D, h, m, s = Tempo.jd2calhms(Tempo.DJ2000, 366.0)
    @test ((s/60 + m)/60 + h)/24 ≈ 0.5

    # output types 
    calhmsec = Tempo.jd2calhms(jd1, jd2+fd-0.5)
    @test all(typeof.(calhmsec) .== (Int, Int, Int, Int, Int, Float64))
    @test all(typeof.((jd1, jd2)) .== Int64)
end

@testset "Function leapseconds" begin 
    for (y, m, l) in Tempo.LEAP_TABLE[1:end]
        @test Tempo.leapseconds(y, m) == l 
        @test Tempo.leapseconds(y, m+1) == l 
        @test Tempo.leapseconds(y, m-1) != l 
    end
end