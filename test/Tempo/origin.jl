@testset "Origin offsets" begin
    iy = 1993
    im = 2
    id = 25
    ihour = 12
    imin = 53
    sec = 15.0

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
