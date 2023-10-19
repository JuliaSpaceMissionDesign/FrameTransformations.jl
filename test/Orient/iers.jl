@testset "IERS Transformations" verbose = true begin
    atol, rtol = 1e-11, 1e-11

    # Radians to arcseconds
    r2a = 180 / π * 3600

    # Function to compute the angle between 2 vectors in arcseconds
    v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600

    @testset "Fundamental Arguments" verbose = true begin

        # Testing Delaunay's Arguments 
        @testset "Fundamental Arguments 2003" begin
            t = [rand(29)..., 1.0]

            for i in eachindex(t)
                fa = Orient.FundamentalArguments(t[i], rand([iau2000a, iau2006a]))

                # --- Delaunay Arguments 
                # Testing Mean Anomaly of the Moon 
                @test fa.Mₐ ≈ fal03(t[i]) atol = atol rtol = rtol

                # Testing Mean Anomaly of the Sun 
                @test fa.Sₐ ≈ falp03(t[i]) atol = atol rtol = rtol

                # Testing Mean Argument of Latitude of the Moon 
                @test fa.uₘ ≈ faf03(t[i]) atol = atol rtol = rtol

                # Testing mean elongation of the moon from the sun 
                @test fa.Dₛ ≈ fad03(t[i]) atol = atol rtol = rtol

                # Testing Mean Longitude of the Moon
                @test fa.Ωₘ ≈ mod(faom03(t[i]), 2π) atol = atol rtol = rtol

                # --- Planetary Arguments 
                # Testing Mean Longitude of Mercury
                @test fa.λ_Me ≈ fame03(t[i]) atol = atol rtol = rtol

                # Testing Mean Longitude of Venus
                @test fa.λ_Ve ≈ fave03(t[i]) atol = atol rtol = rtol

                # Testing Mean Longitude of Earth
                @test fa.λ_Ea ≈ fae03(t[i]) atol = atol rtol = rtol

                # Testing Mean Longitude of Mars
                @test fa.λ_Ma ≈ fama03(t[i]) atol = atol rtol = rtol

                # Testing Mean Longitude of Jupiter
                @test fa.λ_Ju ≈ faju03(t[i]) atol = atol rtol = rtol

                # Testing Mean Longitude of Saturn
                @test fa.λ_Sa ≈ fasa03(t[i]) atol = atol rtol = rtol

                # Testing Mean Longitude of Uranus
                @test fa.λ_Ur ≈ faur03(t[i]) atol = atol rtol = rtol

                # Testing Mean Longitude of Neptune
                @test fa.λ_Ne ≈ fane03(t[i]) atol = atol rtol = rtol

                # Testing General Accumulated Precession 
                @test fa.pₐ ≈ fapa03(t[i]) atol = atol rtol = rtol
            end
        end

        t = [0.06567, 1.0]

        # Testing the truncated expressions of the Delunary Arguments 
        fab = [
            Orient.FundamentalArguments(t[1], rand([iau2000b, iau2006b])),
            Orient.FundamentalArguments(t[2], rand([iau2000b, iau2006b])),
        ]

        fa = [Orient.FundamentalArguments(t[1]), Orient.FundamentalArguments(t[2])]

        @testset "Delaunay Arguments IAU2000B" begin

            # -- Testing Mean Anomaly of the Moon 
            @test fab[1].Mₐ ≈ 2.663599945842266 atol = atol rtol = rtol
            @test fab[2].Mₐ ≈ 5.826449449627781 atol = atol rtol = rtol

            # -- Testing Mean Anomaly of the Sun 
            @test fab[1].Sₐ ≈ 3.518352372771510 atol = atol rtol = rtol
            @test fab[2].Sₐ ≈ 6.223484580361229 atol = atol rtol = rtol

            # -- Testing Mean Argument of Latitude of the Moon 
            @test fab[1].uₘ ≈ 2.533320574432192 atol = atol rtol = rtol
            @test fab[2].uₘ ≈ 3.059379762905646 atol = atol rtol = rtol

            # -- Testing mean elongation of the moon from the sun 
            @test fab[1].Dₛ ≈ 3.236085510636535e-1 atol = atol rtol = rtol
            @test fab[2].Dₛ ≈ 4.275387201216140 atol = atol rtol = rtol

            # -- Testing Mean Longitude of the Moon
            @test fab[1].Ωₘ ≈ -3.438601115926876e-2 + 2π atol = atol rtol = rtol
            @test fab[2].Ωₘ ≈ -1.586802211172697e-1 + 2π atol = atol rtol = rtol

            for f in [:λ_Me, :λ_Ve, :λ_Ea, :λ_Ma, :λ_Ju, :λ_Sa, :λ_Ur, :λ_Ne, :pₐ]
                for i in eachindex(fab)
                    @test getproperty(fab[i], f) == getproperty(fa[i], f)
                end
            end
        end
    end

    @testset "Obliquity" verbose = true begin
        atol, rtol = 1e-8, 1e-8

        t = [rand(0.0:20000, 49)..., 0.0]

        for i in eachindex(t)
            ep = Epoch("$(t[i]) TT")

            tt_d = Tempo.j2000(ep)
            tt_c = Tempo.j2000c(ep)

            ϵ80 = obl80(DJ2000, tt_d) * r2a
            ϵ06 = obl06(DJ2000, tt_d) * r2a

            @test ϵ80 ≈ r2a * orient_obliquity(iau1980, tt_c) atol = atol rtol = rtol
            @test ϵ06 ≈ r2a * orient_obliquity(iau2006a, tt_c) atol = atol rtol = rtol
        end
    end

    @testset "Nutation" verbose = true begin

        # Testing IAU 2000A Nutation model from ERFA 
        # must be precise up to 1 μas
        atol, rtol = 1e-7, 1e-7

        for _ in 1:50
            ep = Epoch("$(rand(0.0:20000)) TT")

            tt_d = Tempo.j2000(ep)
            tt_c = Tempo.j2000c(ep)

            fa = Orient.FundamentalArguments(tt_c)

            # -- Testing IAU2000A Nutation model
            Δψ, Δϵ = Orient.nutation00(iau2000a, tt_c, fa) .* r2a
            p, e = nut00a(DJ2000, tt_d) .* r2a

            @test Δψ ≈ p atol = atol rtol = rtol
            @test Δϵ ≈ e atol = atol rtol = rtol

            # -- Testing IAU2006A Nutation model 
            Δψ, Δϵ = Orient.orient_nutation(iau2006a, tt_c) .* r2a
            p, e = nut06a(DJ2000, tt_d) .* r2a

            @test Δψ ≈ p atol = atol rtol = rtol
            @test Δϵ ≈ e atol = atol rtol = rtol

            # -- Testing IAU 2000B Nutation model from ERFA 
            # Must be precise up to n 1 mas 
            m = rand([iau2000b, iau2006b])
            Δψ, Δϵ = Orient.orient_nutation(m, tt_c) .* r2a
            p, e = nut00b(DJ2000, tt_d) .* r2a

            @test Δψ ≈ p atol = atol rtol = rtol
            @test Δϵ ≈ e atol = atol rtol = rtol
        end
    end

    @testset "Precession" verbose = true begin
        atol, rtol = 1e-9, 1e-9
        t = [rand(0.0:20000, 49)..., 0.0]

        # -- Testing Frame Bias IAU 2000 (does not depend on t)
        Δψb, Δϵb, Δα₀ = Orient.frame_bias(iau2000a) .* r2a
        db, de, da = bi00() .* r2a

        @test Δψb ≈ db atol = atol
        @test Δϵb ≈ de atol = atol
        @test Δα₀ ≈ da atol = atol

        for i in eachindex(t)
            ep = Epoch("$(t[i]) TT")
            tt_d = Tempo.j2000(ep)
            tt_c = Tempo.j2000c(ep)

            m2000 = rand([iau2000a, iau2000b])
            m2006 = rand([iau2006a, iau2006b])

            # -- Testing Fukushima-Williams angles IAU 2006A
            fw = Orient.fw_angles(m2006, tt_c)
            fe = pfw06(DJ2000, tt_d)

            for j in 1:4
                @test fw[j] .* r2a ≈ fe[j] .* r2a atol = atol rtol = rtol
            end

            v = rand(BigFloat, 3)
            v /= norm(v)

            # -- Testing Rotation FW Rotation Matrix IAU 2006A
            Rₑ = fw2m(fe[1], fe[2], fe[3], fe[4])
            Rₐ = Orient.fw_matrix(fw[1], fw[2], fw[3], fw[4])

            @test v2as(Rₑ * v, Rₐ * v) ≤ 1e-10

            # -- Testing Precession Rate IAU 2000
            Δψₚ, Δϵₚ = Orient.precession_rate(m2000, tt_c) .* r2a
            dp, de = pr00(DJ2000, tt_d) .* r2a

            @test Δψₚ ≈ dp atol = atol rtol = rtol
            @test Δϵₚ ≈ de atol = atol rtol = rtol

            # -- Testing Bias-Precession matrix 
            # IAU 2000 
            Ra = Orient.orient_bias_precession(m2000, tt_c)
            Re = pmat00(DJ2000, tt_d)

            @test v2as(Ra * v, Re * v) ≤ 1e-10

            # IAU 2006 
            Ra = Orient.orient_bias_precession(m2006, tt_c)
            Re = pmat06(DJ2000, tt_d)

            @test v2as(Ra * v, Re * v) ≤ 1e-10

            # -- Testing Bias-Precession-Nutation
            # IAU 2000A 
            Ra = Orient.orient_bias_precession_nutation(iau2000a, tt_c)
            Re = pnm00a(DJ2000, tt_d)
            @test v2as(Ra * v, Re * v) ≤ 1e-7

            # IAU 2000B 
            Ra = Orient.orient_bias_precession_nutation(iau2000b, tt_c)
            Re = pnm00b(DJ2000, tt_d)
            @test v2as(Ra * v, Re * v) ≤ 1e-7

            # IAU 2006A 
            Ra = Orient.orient_bias_precession_nutation(iau2006a, tt_c)
            Re = pnm06a(DJ2000, tt_d)
            @test v2as(Ra * v, Re * v) ≤ 1e-7
        end
    end

    @testset "ITRF to GCRF Routines" verbose = true begin
        atol, rtol = 1e-12, 1e-12
        t = [rand(0.0:20000, 49)..., 0.0]

        @testset "Polar Motion" begin
            for i in eachindex(t)
                ep = Epoch("$(t[i]) TT")
                tt_c = Tempo.j2000c(ep)

                # -- Testing TIO Locator 
                sp = Orient.tio_locator(tt_c) .* r2a
                spₑ = sp00(DJ2000, Tempo.j2000(ep)) .* r2a

                @test sp ≈ spₑ atol = atol rtol = rtol

                # -- Testing Polar Motion 
                xₚ, yₚ = rand(), rand(), rand()

                RPₐ = Orient.polar_motion(tt_c, xₚ, yₚ)
                RPₑ = pom00(xₚ, yₚ, sp00(DJ2000, Tempo.j2000(ep)))'

                v = rand(BigFloat, 3)
                v /= norm(v)
                @test v2as(RPₑ * v, RPₐ * v) ≤ 1e-9
            end
        end

        @testset "Earth Rotation Angle" begin
            for i in eachindex(t)
                ep = Epoch("$(t[i]) TT")
                ep_ut1 = convert(UT1, ep)

                ut1_d = Tempo.j2000(ep_ut1)

                # -- Testing ERA Rotation Angle 
                ERA = Orient.earth_rotation_angle(ut1_d) .* r2a
                ERAₑ = era00(DJ2000, ut1_d) .* r2a

                @test ERA ≈ ERAₑ atol = 1e-9 rtol = 1e-9

                # -- Testing Earth velocity 
                @test Orient.earth_rotation_rate() == Orient.earth_rotation_rate(0)
            end
        end

        @testset "Precession-Nutation" begin
            for i in eachindex(t)
                ep = Epoch("$(t[i]) TT")
                tt_d = Tempo.j2000(ep)
                tt_c = Tempo.j2000c(ep)

                v = rand(BigFloat, 3)
                v /= norm(v)

                m2000 = rand([iau2000a, iau2000b])
                m2006 = rand([iau2006a, iau2006b])

                # -- Testing CIP coordinates 
                # IAU 2000A 
                X, Y = Orient.cip_coords(iau2000a, tt_c) .* r2a
                Xₑ, Yₑ = bpn2xy(pnm00a(DJ2000, tt_d)) .* r2a

                @test X ≈ Xₑ atol = 1e-7
                @test Y ≈ Yₑ atol = 1e-7

                # IAU 2000B 
                X, Y = Orient.cip_coords(iau2000b, tt_c) .* r2a
                Xₑ, Yₑ = bpn2xy(pnm00b(DJ2000, tt_d)) .* r2a

                @test X ≈ Xₑ atol = 1e-7
                @test Y ≈ Yₑ atol = 1e-7

                # IAU 2006A 
                X, Y = Orient.cip_coords(iau2006a, tt_c) .* r2a
                Xₑ, Yₑ = bpn2xy(pnm06a(DJ2000, tt_d)) .* r2a

                @test X ≈ Xₑ atol = 1e-7
                @test Y ≈ Yₑ atol = 1e-7

                # -- Testing CIO Locator 
                X, Y = rand(), rand()

                # IAU 2000A/B 
                s = Orient.cio_locator(m2000, tt_c, X, Y) * r2a
                sₑ = s00(DJ2000, tt_d, X, Y) * r2a
                @test s ≈ sₑ atol = 1e-12

                # IAU 2006A/B 
                s = Orient.cio_locator(m2006, tt_c, X, Y) * r2a
                sₑ = s06(DJ2000, tt_d, X, Y) * r2a
                @test s ≈ sₑ atol = 1e-12

                # -- Testing CIP Rotation Matrix 
                # IAU 2000A 
                Q = Orient.cip_motion(iau2000a, tt_c, 0.0, 0.0)
                Qₑ = c2i00a(DJ2000, tt_d)'
                @test v2as(Q * v, Qₑ * v) ≤ 1e-7

                # IAU 2000B 
                Q = Orient.cip_motion(iau2000b, tt_c, 0.0, 0.0)
                Qₑ = c2i00b(DJ2000, tt_d)'
                @test v2as(Q * v, Qₑ * v) ≤ 1e-7

                # IAU 2006A 
                Q = Orient.cip_motion(iau2006a, tt_c, 0.0, 0.0)
                Qₑ = c2i06a(DJ2000, tt_d)'
                @test v2as(Q * v, Qₑ * v) ≤ 1e-7
            end
        end

        @testset "Full Rotation" verbose = true begin
            for i in eachindex(t)
                ep = Epoch("$(t[i]) TT")
                ep_utc = convert(UTC, ep)
                ep_ut1 = convert(UT1, ep_utc)

                tt_d = Tempo.j2000(ep)
                tt_s = Tempo.j2000s(ep)
                tt_c = Tempo.j2000c(ep)

                utc_d = Tempo.j2000(ep_utc)

                # We pass through this interpolation otherwise some tests fail!
                offset = interpolate(Orient.IERS_EOP.UT1_TT, tt_d)
                ut1_d = tt_d + offset / Tempo.DAY2SEC

                v = rand(BigFloat, 3)
                v /= norm(v)

                # Polar coordinates with UTC 
                xₚ = arcsec2rad(interpolate(Orient.IERS_EOP.x_TT, tt_d))
                yₚ = arcsec2rad(interpolate(Orient.IERS_EOP.y_TT, tt_d))

                # Polar coordinates with TT 
                xₚ_TT = arcsec2rad(interpolate(Orient.IERS_EOP.x_TT, tt_d))
                yₚ_TT = arcsec2rad(interpolate(Orient.IERS_EOP.y_TT, tt_d))

                # CIP Deviations 
                dX = arcsec2rad(1e-3 * interpolate(Orient.IERS_EOP.dX_TT, tt_d))
                dY = arcsec2rad(1e-3 * interpolate(Orient.IERS_EOP.dY_TT, tt_d))

                # -- Testing GCRF-to-ITRF Rotation 
                # IAU 2000A (accurate to ~μas)
                Rₑ = c2t00a(DJ2000, tt_d, DJ2000, ut1_d, xₚ, yₚ)'
                R = orient_rot3_itrf_to_gcrf(iau2000a, tt_s, ut1_d, xₚ, yₚ)
                @test v2as(R * v, Rₑ * v) ≤ 1e-7

                # IAU 2000B (accurate to ~mas)
                Rₑ = c2t00b(DJ2000, tt_d, DJ2000, ut1_d, xₚ, yₚ)'
                R = orient_rot3_itrf_to_gcrf(iau2000b, tt_s, ut1_d, xₚ, yₚ)
                @test v2as(R * v, Rₑ * v) ≤ 1e-4

                # IAU 2006A (accurate to ~μas)
                Rₑ = c2t06a(DJ2000, tt_d, DJ2000, ut1_d, xₚ, yₚ)'
                R = orient_rot3_itrf_to_gcrf(iau2006a, tt_s, ut1_d, xₚ, yₚ)
                @test v2as(R * v, Rₑ * v) ≤ 1e-7

                # -- Testing Time Routines 
                # IAU2000A \ IAU2006A
                for m in (iau2000a, iau2006a)
                    R1 = orient_rot3_itrf_to_gcrf(m, tt_s)
                    R2 = orient_rot3_itrf_to_gcrf(m, tt_s, ut1_d, xₚ, yₚ, dX, dY)
                    @test v2as(R1 * v, R2 * v) ≤ 1e-10
                end

                # IAU2000B
                offset = interpolate(Orient.IERS_EOP.UT1_TT, tt_d)
                ut1 = tt_d + offset / Tempo.DAY2SEC

                R1 = orient_rot3_itrf_to_gcrf(iau2000b, tt_s)
                R2 = orient_rot3_itrf_to_gcrf(iau2000b, tt_s, ut1, xₚ_TT, yₚ_TT, 0.0, 0.0)
                @test v2as(R1 * v, R2 * v) ≤ 1e-10

                # CPNc (accurate to ~16mas)
                R = orient_rot3_itrf_to_gcrf(CPNc, tt_s)
                Rₑ = orient_rot3_itrf_to_gcrf(iau2006a, tt_s)
                @test v2as(R * v, Rₑ * v) ≤ 1e-1

                # CPNc (accurate to ~1as)
                R = orient_rot3_itrf_to_gcrf(CPNd, tt_s)
                Rₑ = orient_rot3_itrf_to_gcrf(iau2006a, tt_s)
                @test v2as(R * v, Rₑ * v) ≤ 1
            end
        end
    end
end;
