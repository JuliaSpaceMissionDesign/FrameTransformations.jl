@testset "instant.jl" begin
    s = 1.0seconds
    m = 1.0minutes
    h = 1.0hours
    d = 1.0days
    w = 1.0weeks
    M = 1.0months
    q = 1.0quarters
    y = 1.0years
    c = 1.0centuries

    @test s == Instant(seconds, 1.0)
    @test m == Instant(minutes, 1.0)
    @test h == Instant(hours, 1.0)
    @test d == Instant(days, 1.0)
    @test w == Instant(weeks, 1.0)
    @test M == Instant(months, 1.0)
    @test q == Instant(quarters, 1.0)
    @test y == Instant(years, 1.0)
    @test c == Instant(centuries, 1.0)

    @test seconds(s) == 1.0seconds
    @test seconds(m) == seconds(s) * 60
    @test seconds(h) == seconds(m) * 60
    @test seconds(d) == seconds(h) * 24
    @test seconds(w) == seconds(d) * 7
    @test seconds(M) == seconds(d) * 30
    @test seconds(q) == seconds(M) * 4
    @test seconds(y) == 3.15576e7seconds
    @test seconds(c) == 3.15576e9seconds

    @test zero(Instant{AstronautBase.Tempo.Year}) == 0.0years
    @test zero(1years) == 0years
    @test zero(1.0years) == 0.0years

    @test eltype(1.0years) == Float64
    @test eltype(typeof(1.0years)) == Float64

    a = [1.0, 2.0, 3.0]
    @test a * seconds == [1.0seconds, 2.0seconds, 3.0seconds]
    @test a .* seconds == [1.0seconds, 2.0seconds, 3.0seconds]
    @test seconds * a == [1.0seconds, 2.0seconds, 3.0seconds]
    @test seconds .* a == [1.0seconds, 2.0seconds, 3.0seconds]
    @test Base.broadcastable(seconds) isa typeof(Ref(seconds))

    @test unit(1years) == years
    int_rng = 1seconds:3seconds
    @test step(int_rng) == 1seconds
    @test collect(int_rng) == [1seconds, 2seconds, 3seconds]
    float_rng = 1.0seconds:3.0seconds
    @test step(float_rng) == 1.0seconds
    @test collect(float_rng) == [1.0seconds, 2.0seconds, 3.0seconds]
    @test Instant{Tempo.Second,Float64}(1.0seconds) == 1.0seconds

    @test Tempo.name(seconds) == "seconds"
    @test Tempo.name(minutes) == "minutes"
    @test Tempo.name(hours) == "hours"
    @test Tempo.name(days) == "days"
    @test Tempo.name(weeks) == "weeks"
    @test Tempo.name(months) == "months"
    @test Tempo.name(quarters) == "quarters"
    @test Tempo.name(years) == "years"
    @test Tempo.name(centuries) == "centuries"

end