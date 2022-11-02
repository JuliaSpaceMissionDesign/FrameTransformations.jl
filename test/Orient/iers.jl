
@testset "Fundamental Arguments" begin 

    # Testing Fundamental Arguments 
    # 1st set data from Vallado 4th ed. Example 3-14 pp. 220
    # 2nd set data from ERFA values in radians (guarantees a higher precision!)
    
    # Missing the longitude of neptune as Vallado uses a different value wrt to 
    # SOFA\IERS 2010 conventions 

    tv = 0.0426236319
    te = 0.12345

    fav = FundamentalArguments(tv)
    fae = FundamentalArguments(te)

    # Testing Delaunay's Arguments 
    @testset "Delaunay arguments 2003" begin

        # Testing Mean Anomaly of the Moon 
        @test rad2deg(fav.Mₐ) ≈ 314.9122873 atol=1e-7 rtol=1e-7
        @test fae.Mₐ ≈ 9.0124422693240238e-2 atol=1e-6

        # Testing Mean Anomaly of the Sun 
        @test rad2deg(fav.Sₐ) ≈ 91.9393769 atol=1e-7 rtol=1e-7
        @test fae.Sₐ ≈ 2.122527458615575 atol=1e-9

        # Testing Mean Argument of Latitude of the Moon 
        @test rad2deg(fav.uₘ) ≈ 169.0970043 atol=1e-7 rtol=1e-7
        @test fae.uₘ ≈ 6.013725526098429 atol=1e-9

        # Testing mean elongation of the moon from the sun 
        @test rad2deg(fav.Dₛ) ≈ 196.7516428 atol=1e-7 rtol=1e-7
        @test fae.Dₛ ≈ 3.247622743094825 atol=1e-9

        # Testing Mean Longitude of the Moon
        @test rad2deg(fav.Ωₘ) ≈ 42.6046467 atol=1e-7 rtol=1e-7
        @test fae.Ωₘ ≈ -1.984867574205391+2π atol=1e-9 

    end

    @testset "Delaunay arguments IAU2000B" begin 
        # Testing the truncated expressions of the Delunary Arguments 

        fa = FundamentalArguments(te, iau2006b);

        # Testing Mean Anomaly of the Moon 
        @test fa.Mₐ ≈ 9.012187106285367e-02 atol=1e-6

        # Testing Mean Anomaly of the Sun 
        @test fa.Sₐ ≈ 2.122527499497290 atol=1e-9

        # Testing Mean Argument of Latitude of the Moon 
        @test fa.uₘ ≈ 6.013726468232121 atol=1e-9

        # Testing mean elongation of the moon from the sun 
        @test fa.Dₛ ≈ 3.247623213717811 atol=1e-9

        # Testing Mean Longitude of the Moon
        @test fa.Ωₘ ≈ -1.984868126360060+2π atol=1e-9 

        for f in [:λ_Me, :λ_Ve, :λ_Ea, :λ_Ma, :λ_Ju, :λ_Sa, :λ_Ur, :λ_Ne, :pₐ]
            @test getproperty(fa, f) == getproperty(fae, f)
        end

    end 

    # Testing Planetary Arguments 
    @testset "Planetary Arguments" begin

        # Testing Mean Longitude of Mercury
        @test rad2deg(fav.λ_Me) ≈ 143.319167 atol=1e-7 rtol=1e-7
        @test fae.λ_Me ≈ 6.015322458572193 atol=1e-14 rtol=1e-14

        # Testing Mean Longitude of Venus
        @test rad2deg(fav.λ_Ve) ≈ 156.221635 atol=1e-7 rtol=1e-7
        @test fae.λ_Ve ≈ 3.595450621383065 atol=1e-14 rtol=1e-14

        # Testing Mean Longitude of Earth
        @test rad2deg(fav.λ_Ea) ≈ 194.890465 atol=1e-7 rtol=1e-7
        @test fae.λ_Ea ≈ 3.919817995983848 atol=1e-14 rtol=1e-14

        # Testing Mean Longitude of Mars
        @test rad2deg(fav.λ_Ma) ≈ 91.262347 atol=1e-7 rtol=1e-7
        @test fae.λ_Ma ≈ 3.461044170354398 atol=1e-14 rtol=1e-14

        # Testing Mean Longitude of Jupiter
        @test rad2deg(fav.λ_Ju) ≈ 163.710186 atol=1e-7 rtol=1e-7
        @test fae.λ_Ju ≈ 8.553961236235601e-1 atol=1e-14 rtol=1e-14

        # Testing Mean Longitude of Saturn
        @test rad2deg(fav.λ_Sa) ≈ 102.168400 atol=1e-7 rtol=1e-7
        @test fae.λ_Sa ≈ 3.507194207731200 atol=1e-14 rtol=1e-14

        # Testing Mean Longitude of Uranus
        @test rad2deg(fav.λ_Ur) ≈ 332.317825 atol=1e-7 rtol=1e-7
        @test fae.λ_Ur ≈ 1.212873991300292e-1 atol=1e-14 rtol=1e-14

        # Testing Mean Longitude of Neptune (only from SOFA routines)
        @test fae.λ_Ne ≈ 5.782638611951110 atol=1e-20 rtol=1e-20

        # Testing General Accumulated Precession 
        @test rad2deg(fav.pₐ) ≈ 0.059545 atol=1e-6 rtol=1e-7
        @test fae.pₐ ≈ 3.010009133483176e-3 atol=1e-14 rtol=1e-14

    end 

end;

@testset "nutation" begin 

    t = -0.038913073237503454; 

    ERFA_DJ00 = 2451545.
    ERFA_DJC = 36525.

    date1 = 2450123.7;
    t = ((date1 - ERFA_DJ00)) / ERFA_DJC;

    # Testing IAU 2000B Nutation model from ERFA 

    fa = FundamentalArguments(t, iau2006b);
    Δψ, Δϵ = nutation00(iau2006b, t, fa)

    @test Δψ ≈  3.545257283580251e-05 atol = 1e-8
    @test Δϵ ≈ -4.139189402616098e-05 atol = 1e-8

    # Testing IAU 2000B Nutation model from ERFA 

    fa = FundamentalArguments(t);
    Δψ, Δϵ = nutation00(iau2006a, t, fa)

    @test Δψ ≈  3.545257283580251e-05 atol = 1e-8
    @test Δϵ ≈ -4.139189402616098e-05 atol = 1e-8
    
end

@testset "itrf2gcrf" begin

    # Testing TIO Locator (from Vallado 4th ed. Example 3-14, pp. 220)
    TT = 0.0426236319
    sp = tio_locator(TT)

    @test sp ≈ -9.712e-12 atol=1e-12 rtol=1e-12

    # Testing Polar Motion Matrix (from Vallado 4th ed. Example 3-14, pp. 220)
    xₚ = -0.140682 |> arcsec2rad
    yₚ = 0.333309 |> arcsec2rad
    TT = 0.0426236319
    
    W = polar_motion(xₚ, yₚ, tio_locator(TT))

    r_ITRF = SA[-1033.4793830, 7901.2952754, 6380.3565958]
    r_TIRS = SA[-1033.4750312, 7901.3055856, 6380.3445327]

    @test maximum(abs.(W*r_ITRF - r_TIRS)) ≈ 0. atol=1e-7 rtol=1e-7

    # Testing ERA (from Vallado 4th ed. Example 3-14, pp. 221)
    Tᵤ = 2453101.827406783 - 2451545
    θₑ = deg2rad(312.7552829)

    @test earth_rotation_angle(Tᵤ) ≈ θₑ atol=1e-7 rtol=1e-7

    # Testing ERA Matrix (from Vallado 4th ed. Example 3-14, pp. 221)
    r_TIRS = SA[-1033.4750312, 7901.3055856, 6380.3445327]
    r_CIRS = SA[5100.0184047, 6122.7863648, 6380.3445327]

    Tᵤ = 2453101.827406783 - 2451545
    Rₑ = era_rotm(Tᵤ)

    @test maximum(abs.(Rₑ*r_TIRS - r_CIRS)) ≈ 0. atol=1e-7 rtol=1e-7

    
end

