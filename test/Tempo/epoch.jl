@testset "epoch.jl" begin
    
    dt = DateTime(3600.0)
    e = Epoch(dt, Tempo.TDB)

    @test value(e) ≈ 3600
    @test timescale(e) == Tempo.TDB 

    @test value(e-60minutes) ≈ 0.0
    @test value(e+1hours)seconds ≈ seconds(2hours)
    @test value(e-3600) ≈ 0.0
    @test value(e+3600)seconds ≈ seconds(2hours)

    @test e ≈ Epoch("2000-01-01T13:00:00.0", TDB)
    @test e ≈ Epoch(3600.0, TDB)
    @test e ≈ Epoch{BarycentricDynamicalTime}(3600.0)
    @test e ≈ Epoch{Tempo.BarycentricDynamicalTime}(3600, 0.0)
    @test collect(e:1hours:(e+1hours)) == [e, e+1hours]
end