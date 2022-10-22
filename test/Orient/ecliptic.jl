@testset "ecleod" begin

    # computed with ERFA
    MrefT = SMatrix{3,3}(1.00000000e+00, -7.07836896e-08,  8.05621398e-08,
                         3.28970041e-08,  9.17482130e-01,  3.97776999e-01,
                         -1.02070447e-07, -3.97776999e-01,  9.17482130e-01)
    M = orient_icrs2ecleod(2400000.5, 51544.0)
    @test all(isapprox.(M*MrefT - I(3), 0.0, atol=1e-6))

end

