
@axes ICRF 1 InternationalCelestialReferenceFrame
@axes IAU_TEST 2
@axes IAU_MIMAS 10039

@point MIMAS 601

# load constants 
kclear()

@testset "Body-Centered Rotating (TOD)" verbose=true begin
    
    function mimas(sec)

        D = sec/Tempo.DAY2SEC
        T = D/Tempo.CENTURY2DAY

        S3 = 177.40 - 36505.5*T
        S5 = 316.45 + 506.2*T

        a = 40.66 - 0.036*T + 13.56*sind(S3)
        d = 83.52 - 0.004*T - 1.53*cosd(S3)
        w = 333.46 + 381.994555*D  - 13.48*sind(S3) - 44.85*sind(S5)

        return deg2rad.([a, d, w])
    end

    function R_mimas(sec)
        ra, dec, w = mimas(sec)
        angle_to_dcm(π/2 + ra, π/2 - dec, w, :ZXZ)
    end 

    tpc10_constants = Basic.load(TPC(path(KERNELS[:IAU10])));
    furnsh(path(KERNELS[:IAU10]));
    
    # Test against manual computation 
    p = Basic.Orient.PlanetsPrecessionNutation(point_alias(MIMAS), tpc10_constants)
    f1, f2, f3 = Basic.Orient.orient_planets_angles(p, "mimas")


    @testset "IAU Euler angles" begin
        
        for _ in 1:50
            e = rand(0.0:1e8)

            # angle 
            angles = [f1(e)...]
            angles_true = mimas(e)
            @test angles ≈ angles_true atol=1e-14 rtol=1e-14

            # 1st derivative 
            δangles = [f2(e)[4:end]...]
            δangles_true = D¹(mimas, e)
            @test δangles ≈ δangles_true atol=1e-14 rtol=1e-14

            # 2nd derivative 
            δ²angles = [f3(e)[7:end]...]
            δ²angles_true = D²(mimas, e)
            @test δ²angles ≈ δ²angles_true atol=1e-14 rtol=1e-14
        end 

    end;

    @testset "2nd/3rd order AD with MIMAS" verbose=false begin 

        v2as = (x, y) -> acosd(max(-1, min(1, dot(x/norm(x), y/norm(y)))))*3600

        # Create dummy frame system 
        FRAMES = FrameSystem{4, Float64}()

        add_axes_inertial!(FRAMES, ICRF)
        add_axes_bcrtod!(FRAMES, tpc10_constants, MIMAS, IAU_MIMAS, ICRF)

        for _ in 1:50
            e = rand(0.0:1e8)
            Rb = rotation12(FRAMES, ICRF, IAU_MIMAS, e)

            v = rand(BigFloat, 3);
            v /= norm(v)

            # -- Test 2nd derivative
            Re = D²(R_mimas, e)
            @test v2as(Re*v, Rb[3]*v) ≤ 1e-6
            
            # -- Test 3rd derivative
            Re = D³(R_mimas, e)
            @test v2as(Re*v, Rb[4]*v) ≤ 1e-6

        end 
    end;

    for (version, pckname) in zip((10, 11), (:IAU10, :IAU11))
        # load constants 
        kclear()

        tpc_constants = Basic.load(TPC(path(KERNELS[pckname]))); 
        furnsh(path(KERNELS[pckname]))

        NAIFIds = []
        for k in keys(tpc_constants)
            (k>10 && k<1000 && haskey(tpc_constants[k], :pole_ra)) ? push!(NAIFIds, k) : ()
        end

        v2as = (x, y) -> acosd(max(-1, min(1, dot(x/norm(x), y/norm(y)))))*3600

        @testset "1st and 2nd order Rotations PCK$version" verbose=false begin 
            @testset "Body with NAIFId $i" for i in NAIFIds 

                FRAMES = FrameSystem{3, Float64}()
                add_axes_inertial!(FRAMES, ICRF)
                Frames._axes_bcrtod!(FRAMES, tpc_constants, i, "test", IAU_TEST, ICRF)

                for _ in 1:25
                    ep = rand(0.0:1e6)

                    R = sxform("J2000", "IAU_$(bodc2n(i))", ep)
                    R_ = rotation6(FRAMES, ICRF, IAU_TEST, ep)

                    v = rand(BigFloat, 3)
                    v /= norm(v)

                    # It has a precision of about 1 μas
                    @test v2as(R[1:3, 1:3]*v, R_[1]*v) ≤ 1e-6
                    @test v2as(R[4:6, 1:3]*v, R_[2]*v) ≤ 1e-6
                    
                end
            end;
        end;

        # unload kernels!
        kclear()
    end

end
