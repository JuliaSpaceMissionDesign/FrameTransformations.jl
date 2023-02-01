kclear()

@axes ICRF 1 InternationalCelestialReferenceFrame
@axes ECLIPJ2000 17
@axes ITRF 3000

@testset "ITRF" verbose=true begin 

    v2as = (x, y) -> acosd(max(-1, min(1, dot(x/norm(x), y/norm(y)))))*3600

    # Check that if you haven't loaded the kernel you get an error 
    frames = FrameSystem{2, Float64}() 
    add_axes_inertial!(frames, ICRF)

    @test_throws ErrorException add_axes_pa421!(frames, PA421, ICRF)

    eph = CalcephProvider(path(KERNELS[:ITRF]))

    frames = FrameSystem{2, Float64}(eph) 
    add_axes_inertial!(frames, ICRF)
    add_axes_eclipj2000!(frames, ECLIPJ2000, ICRF)
    add_axes_ephemeris!(frames, ITRF, :ZXZ)

    furnsh(path(KERNELS[:LEAP]), path(KERNELS[:ITRF]))

    for _ = 1:25 
        # FIXME: la data non la considera valida!
        et = 2.48e6 - DJ2000

        v = rand(BigFloat, 3)
        v /= norm(v) 

        Rs = sxform("J2000", "ITRF93", et)
        Rb = rotation6(frames, ICRF, ITRF, et)

        @test v2as(Rs[1:3, 1:3]*v, Rb[1]*v) ≤ 1e-6

    end


end;

kclear()