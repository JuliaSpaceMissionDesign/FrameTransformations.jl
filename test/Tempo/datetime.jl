@testset "datetime.jl" begin
    d = Date(2000, 1, 1)

    @test year(d) == 2000
    @test month(d) == 1 
    @test day(d) == 1
    @test calendar(d) == :gregorian
    @test isleapyear(d) 
    @test finddayinyear(d) == 1 

    @test d == Date(0)
    @test d == Date(2000, 1)
    @test j2000(d) == 0

    @test Tempo.findmonth(32, true) == 2
    @test Tempo.findday(32, 2, true) == 1
    @test d + 1days == Date(1)
    @test d + 86400seconds == Date(1)
    @test Date(2) - 2days == d 
    @test d + 1 == Date(1)
    @test Date(1) - 1 == d 

    t0 = Time(12, 0, 0.0)
    @test hour(t0) == 12 
    @test minute(t0) == 0 
    @test second(t0) == 0
    @test findfractionofday(t0) == 0.5 
    @test findsecondinday(t0) == SECONDS_PER_DAY/2

    t = Time(12, 1, 15.300300300)
    @test millisecond(t) == 300
    @test microsecond(t) == 300
    @test nanosecond(t) == 300
    @test second(Float64, t) == 15.300300300
    @test second(Int64, t) == 15
    @test findfractionofsecond(t) ≈ 0.300300300

    J2000 = DateTime(d, SECONDS_PER_DAY/2)
    @test t0 == Time(J2000)
    @test d == Date(J2000)

    D = DateTime(d, 0.0)
    @test J2000 - D == 43200seconds
    @test J2000 - 43200 == D
    @test J2000 - 43200seconds == D
    @test J2000 - 0.5days == D
    @test J2000 + 0.5days == D + 1days

    s = "2000-01-01T12:00:00.0"
    @test J2000 == DateTime(s)
    @test J2000 == DateTime(0.0)
    
    D2 = DateTime(SECONDS_PER_DAY/3)
    @test j2000(D2)-0.5 ≈ 1/3
    @test j2000seconds(D2) ≈ SECONDS_PER_DAY/3 + SECONDS_PER_DAY/2

    DTR = J2000:(J2000+2days)
    @test DTR == [J2000, J2000+1days, J2000+2days]
end