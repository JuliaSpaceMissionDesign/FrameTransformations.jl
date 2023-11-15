
@testset "GMST" verbose=true begin 

    # Radians to arcseconds
    r2a = 180 / π * 3600

    @testset "GMST" verbose=true begin 

        for _ in 1:50 

            ep = Epoch("$(rand(0.0:20000))")
            ep_ut1 = convert(UT1, ep)

            tt_d = Tempo.j2000(ep)
            tt_c = Tempo.j2000c(ep)

            ut1_d = Tempo.j2000(ep_ut1)
            ut1_c = Tempo.j2000c(ep_ut1)

            # Compute ERA 
            θ = Orient.earth_rotation_angle(ut1_d)

            # Test GMST with IAU2006 models
            b = gmst06(DJ2000, ut1_d, DJ2000, tt_d)
            for m in (iau2006a, iau2006b)
                a = Orient.gmst(m, tt_c, θ)
                @test a ≈ b atol=1e-12 rtol=1e-12
            end 

            # Test GMST with IAU2000 models
            b = gmst00(DJ2000, ut1_d, DJ2000, tt_d)
            for m in (iau2000a, iau2000b)
                a = Orient.gmst(m, tt_c, θ)
                @test a ≈ b atol=1e-12 rtol=1e-12
            end 

            # Test GMST with IAU1980 model 
            # TODO: small loss of accuracy here, why?
            a = Orient.gmst(iau1980, ut1_c)
            b = gmst82(DJ2000, ut1_d)
            @test a ≈ b atol=1e-10 rtol=1e-10

        end


    end

    @testset "Equinoxes Equation" verbose=true begin 

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

            @test a ≈ b atol=1e-7 rtol=1e-7

            # -- Testing IAU2000B Equation of the Equinoxes 
            a = Orient.equinoxes_equation(iau2000b, tt_c)*r2a
            b = ee00b(DJ2000, tt_d)*r2a

            @test a ≈ b atol=1e-7 rtol=1e-7

        end

    end

end;