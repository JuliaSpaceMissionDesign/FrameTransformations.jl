kclear()

@testset "Moon Orientation" verbose=true begin
    
    v2as = (x, y) -> acosd(max(-1, min(1, dot(x/norm(x), y/norm(y)))))*3600

    v = rand(BigFloat, 3)
    v /= norm(v) 

    @testset "DE421" verbose=false begin 

        for kernel in (:LEAP, :PA421, :FK_DE421)
            furnsh(path(KERNELS[kernel]))
        end

        ephempty = CalcephProvider(path(KERNELS[:DE432]))
        @test_throws ErrorException Orient.orient_rot3_icrf_to_pa421(ephempty, 0.0)

        eph = CalcephProvider(path(KERNELS[:PA421]))

        for _ = 1:10
            et = rand(0.0:1e8)

            # Test orientation between ME421 and PA421
            Rs = pxform("MOON_PA", "MOON_ME", et)
            @test v2as(Rs*v, Orient.DCM_MOON_PA421_TO_ME421*v) ≤ 1e-6

            # Test orientation between PA421 and ICRF 
            Rs = pxform("J2000", "MOON_PA", et)
            Rb = Orient.orient_rot3_icrf_to_pa421(eph, et)

            @test v2as(Rs*v, Rb*v) ≤ 1e-6

            # Test with epoch! 
            ep = Epoch("$(rand(1:1e4)) TDB")

            R1 = Orient.orient_rot3_icrf_to_pa421(eph, ep)
            R2 = pxform("J2000", "MOON_PA", j2000s(ep))
            @test v2as(R1*v, R2*v) ≤ 1e-6

        end
    end

    @testset "DE440" verbose=false begin 

        for kernel in (:LEAP, :PA440, :FK_DE440)
            furnsh(path(KERNELS[kernel]))
        end

        ephempty = CalcephProvider(path(KERNELS[:DE432]))
        @test_throws ErrorException Orient.orient_rot3_icrf_to_pa440(ephempty, 0.0)

        eph = CalcephProvider(path(KERNELS[:PA440]))

        for _ = 1:10
            et = rand(0.0:1e8)

            # Test orientation between ME421 and PA440
            Rs = pxform("MOON_PA", "MOON_ME", et)
            @test v2as(Rs*v, Orient.DCM_MOON_PA440_TO_ME421*v) ≤ 1e-6

            # Test orientation between PA440 and ICRF 
            Rs = pxform("J2000", "MOON_PA", et)
            Rb = Orient.orient_rot3_icrf_to_pa440(eph, et)

            @test v2as(Rs*v, Rb*v) ≤ 1e-6

            # Test with epoch! 
            ep = Epoch("$(rand(1:1e4)) TDB")

            R1 = Orient.orient_rot3_icrf_to_pa440(eph, ep)
            R2 = pxform("J2000", "MOON_PA", j2000s(ep))
            @test v2as(R1*v, R2*v) ≤ 1e-6

        end
    end

end;

kclear()