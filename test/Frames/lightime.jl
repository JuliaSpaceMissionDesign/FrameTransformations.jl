# unload spice kernels 
kclear()

# Register axes
@axes ICRF 1 InternationalCelestialReferenceFrame
@axes ECLIPJ2000 17
@axes IAU_EARTH 10013

# Register points 
@point SSB 0 SolarSystemBarycenter
@point EMB 3
@point Sun 10
@point Earth 399
@point SaturnB 6 SaturnBarycenter

# -- Position LT corrections 
@testset "Light Time Corrections" verbose = true begin

    # load spice kernels 
    furnsh(path(KERNELS[:DE432]), path(KERNELS[:PCK10]), path(KERNELS[:LEAP]))

    # Spice ECLIPJ2000 has an older orientation!
    DCM_ECLIPJ2000 = DCM(pxform("J2000", "ECLIPJ2000", 0.0))

    # Load TPC 
    tpc_constants = FrameTransformations.load(TPC(path(KERNELS[:PCK10])))

    # Create frame system
    eph = EphemerisProvider(path(KERNELS[:DE432]))
    FRAMES = FrameSystem{3,Float64}(eph)

    # add axes
    add_axes_inertial!(FRAMES, ICRF)
    add_axes_fixedoffset!(FRAMES, ECLIPJ2000, ICRF, DCM_ECLIPJ2000)
    add_axes_bcrtod!(FRAMES, IAU_EARTH, Earth, tpc_constants)

    # add points
    add_point_root!(FRAMES, SSB, ICRF)
    add_point_ephemeris!(FRAMES, Sun)
    add_point_ephemeris!(FRAMES, EMB)
    add_point_ephemeris!(FRAMES, Earth)
    add_point_ephemeris!(FRAMES, SaturnB)

    tol = 1e-11

    epb = Epoch("2022-01-01T12:00:10.0000 TDB")
    eps = str2et("1 January 2022 12:00:10.0000 TDB")

    # --- LT corrections with 1 iteration
    @testset "LT Corrections" verbose = true begin
        for (bfun, sfun) in zip((vector3, vector6), (spkpos, spkezr))
            for (bcorr, scorr) in zip((LightTime, PlanetaryAberration), ("LT", "LT+S"))
                for (bdir, sdir) in zip((-1, 1), ("", "X"))

                    # Corrections with epoch
                    x = bfun(FRAMES, Sun, SaturnB, ICRF, epb, bcorr, bdir)
                    s = sfun("6", eps, "J2000", sdir * scorr, "SUN")[1]
                    @test x ≈ s atol = tol rtol = tol

                    for _ in 1:25
                        et = rand(0.0:1e8)

                        # Corrections without rotations
                        x = bfun(FRAMES, Sun, SaturnB, ICRF, et, bcorr, bdir)
                        s = sfun("6", et, "J2000", sdir * scorr, "SUN")[1]
                        @test x ≈ s atol = tol rtol = tol

                        # Corrections with inertial rotations
                        x = bfun(FRAMES, Sun, SaturnB, ECLIPJ2000, et, bcorr, bdir)
                        s = sfun("6", et, "ECLIPJ2000", sdir * scorr, "SUN")[1]
                        @test x ≈ s atol = tol rtol = tol

                        # # Corrections with time rotations (AXESCENTER != FROM != TO)
                        x = bfun(
                            FRAMES,
                            Sun,
                            SaturnB,
                            IAU_EARTH,
                            et,
                            bcorr,
                            bdir;
                            axescenter=Earth,
                        )
                        s = sfun("6", et, "IAU_EARTH", sdir * scorr, "SUN")[1]
                        @test x ≈ s atol = tol rtol = tol

                        # # Corrections with time rotations (AXESCENTER = TO)
                        x = bfun(
                            FRAMES, Sun, Earth, IAU_EARTH, et, bcorr, bdir; axescenter=Earth
                        )
                        s = sfun("EARTH", et, "IAU_EARTH", sdir * scorr, "SUN")[1]
                        @test x ≈ s atol = tol rtol = tol

                        # # Corrections with time rotations (AXESCENTER = FROM)
                        x = bfun(
                            FRAMES,
                            Earth,
                            SaturnB,
                            IAU_EARTH,
                            et,
                            bcorr,
                            bdir;
                            axescenter=Earth,
                        )
                        s = sfun("6", et, "IAU_EARTH", sdir * scorr, "EARTH")[1]
                        @test x ≈ s atol = tol rtol = tol
                    end
                end
            end
        end
    end

    # --- LT corrections with 3 iterations
    @testset "LT Converged Corrections" verbose = true begin
        for (bfun, sfun) in zip((vector3, vector6), (spkpos, spkezr))
            for (bcorr, scorr) in zip((LightTime, PlanetaryAberration), ("CN", "CN+S"))
                for (bdir, sdir) in zip((-1, 1), ("", "X"))

                    # Corrections with epoch
                    x = bfun(FRAMES, Sun, SaturnB, ICRF, epb, bcorr, bdir; iters=3)
                    s = sfun("6", eps, "J2000", sdir * scorr, "SUN")[1]
                    @test x ≈ s atol = tol rtol = tol

                    for _ in 1:25
                        et = rand(0.0:1e8)

                        # Corrections without rotations
                        x = bfun(FRAMES, Sun, SaturnB, ICRF, et, bcorr, bdir; iters=3)
                        s = sfun("6", et, "J2000", sdir * scorr, "SUN")[1]
                        @test x ≈ s atol = tol rtol = tol

                        # Corrections with inertial rotations
                        x = bfun(FRAMES, Sun, SaturnB, ECLIPJ2000, et, bcorr, bdir; iters=3)
                        s = sfun("6", et, "ECLIPJ2000", sdir * scorr, "SUN")[1]
                        @test x ≈ s atol = tol rtol = tol

                        # # Corrections with time rotations (AXESCENTER != FROM != TO)
                        x = bfun(
                            FRAMES,
                            Sun,
                            SaturnB,
                            IAU_EARTH,
                            et,
                            bcorr,
                            bdir;
                            axescenter=Earth,
                            iters=3,
                        )

                        s = sfun("6", et, "IAU_EARTH", sdir * scorr, "SUN")[1]
                        @test x ≈ s atol = tol rtol = tol

                        # # Corrections with time rotations (AXESCENTER = TO)
                        x = bfun(
                            FRAMES,
                            Sun,
                            Earth,
                            IAU_EARTH,
                            et,
                            bcorr,
                            bdir;
                            axescenter=Earth,
                            iters=3,
                        )

                        s = sfun("EARTH", et, "IAU_EARTH", sdir * scorr, "SUN")[1]
                        @test x ≈ s atol = tol rtol = tol

                        # # Corrections with time rotations (AXESCENTER = FROM)
                        x = bfun(
                            FRAMES,
                            Earth,
                            SaturnB,
                            IAU_EARTH,
                            et,
                            bcorr,
                            bdir;
                            axescenter=Earth,
                            iters=3,
                        )

                        s = sfun("6", et, "IAU_EARTH", sdir * scorr, "EARTH")[1]
                        @test x ≈ s atol = tol rtol = tol
                    end
                end
            end
        end
    end
end;

# unload spice kernels 
kclear()
