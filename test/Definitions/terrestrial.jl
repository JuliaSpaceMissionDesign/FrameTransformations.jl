
@axes ICRF 1 InternationalCelestialReferenceFrame
@axes ECI -200
@axes ITRF 10 
@axes ITRF2 100

@testset "ITRF" verbose=true begin 
    v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600

    frames = FrameSystem{3,Float64}()
    add_axes_inertial!(frames, ECI)

    # Check that the parent must be the allowed one
    @test_throws ArgumentError add_axes_itrf!(frames, ITRF, ECI)

    frames = FrameSystem{2,Float64}()
    add_axes_inertial!(frames, ICRF)

    add_axes_itrf!(frames, ITRF, ICRF)
    add_axes_itrf!(frames, ITRF2, ICRF, CPNd)

    for j = 1:10
        ep = rand(0.0:1e7)
        for (model, axes) in zip((iers2010b, CPNd), (ITRF, ITRF2))

            R = rotation6(frames, ICRF, axes, ep)
            R_ = Rotation{2}(iers_rot6_gcrf_to_itrf(ep, model)...)

            v = rand(BigFloat, 3)
            v /= norm(v)

            # Test position 
            @test v2as(R[1]*v, R_[1]*v) ≈ 0.0 atol = 1e-14 rtol = 1e-14
            
            # Test velocity
            b = rand(BigFloat, 6)
            RRinv = R*inv(R_)
            @test det(RRinv[1]) .≈ 1.0 atol=1e-14 rtol=1e-14
            @test det(RRinv[2]) .≈ 0.0 atol=1e-14 rtol=1e-14
            @test maximum(abs.((R*b - R_*b))) .≈ 0.0 atol=1e-14 rtol=1e-14
        end
    end

end;