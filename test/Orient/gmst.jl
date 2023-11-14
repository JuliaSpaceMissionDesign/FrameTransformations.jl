
@testset "GMST" verbose=true begin 

    # Radians to arcseconds
    r2a = 180 / π * 3600

    @testset "Equinoxes Equation" verbose = true begin 

        # Testing  Equinoxes equations models from ERFA, they are 
        # precise up to ~ 1μas 
        for _ in 1:50
            ep = Epoch("$(rand(0.0:20000)) TT")

            tt_d = Tempo.j2000(ep)
            tt_c = Tempo.j2000c(ep)

            fa = Orient.FundamentalArguments(tt_c)

            # -- Testing IAU2000 EE Complementary terms 
            a = Orient.ee_complementary(iau2000a, tt_c, fa)*r2a
            b = eect00(DJ2000, tt_d)*r2a

            @test a ≈ b atol=1e-11 rtol=1e-11

            # -- Testing IAU2000A Equation of the Equinoxes 
            a = Orient.equinoxes_equation(iau2000a, tt_c)*r2a
            b = ee00a(DJ2000, tt_d)*r2a

            @test a ≈ b atol=1e-8 rtol=1e-7

            # -- Testing IAU2000B Equation of the Equinoxes 
            a = Orient.equinoxes_equation(iau2000b, tt_c)*r2a
            b = ee00b(DJ2000, tt_d)*r2a

            @test a ≈ b atol=1e-8 rtol=1e-7

        end

    end

end;