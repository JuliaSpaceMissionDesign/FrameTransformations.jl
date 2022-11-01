
@testset "Fundamental Arguments" begin 
    # Testing Fundamental Arguments (from Vallado 4th ed. Example 3-14 pp. 220)
    # Missing the longitude of neptune as Vallado uses a different value wrt to 
    # SOFA\IERS 2010 conventions 

    t = 0.0426236319

    fa = FundamentalArguments(t)

    # Testing Delaunay's Arguments 
    @testset "Delaunay arguments" begin
        # Testing mean elongation of the moon from the sun 
        @test rad2deg(fa.Dₛ) ≈ 196.7516428 atol=1e-7 rtol=1e-7
        # Testing Mean Anomaly of the Moon 
        @test rad2deg(fa.Mₐ) ≈ 314.9122873 atol=1e-7 rtol=1e-7
        # Testing Mean Anomaly of the Sun 
        @test rad2deg(fa.Sₐ) ≈ 91.9393769 atol=1e-7 rtol=1e-7
        # Testing Mean Longitude of the Moon
        @test rad2deg(fa.Ωₘ) ≈ 42.6046467 atol=1e-7 rtol=1e-7
        # Testing Mean Argument of Latitude of the Moon 
        @test rad2deg(fa.uₘ) ≈ 169.0970043 atol=1e-7 rtol=1e-7
    end

    # Testing Planetary Arguments 
    @testset "Planetary Arguments" begin
        # Testing Mean Longitude of Mercury
        @test rad2deg(fa.λ_Me) ≈ 143.319167 atol=1e-7 rtol=1e-7
        # Testing Mean Longitude of Venus
        @test rad2deg(fa.λ_Ve) ≈ 156.221635 atol=1e-7 rtol=1e-7
        # Testing Mean Longitude of Earth
        @test rad2deg(fa.λ_Ea) ≈ 194.890465 atol=1e-7 rtol=1e-7
        # Testing Mean Longitude of Mars
        @test rad2deg(fa.λ_Ma) ≈ 91.262347 atol=1e-7 rtol=1e-7
        # Testing Mean Longitude of Jupiter
        @test rad2deg(fa.λ_Ju) ≈ 163.710186 atol=1e-7 rtol=1e-7
        # Testing Mean Longitude of Saturn
        @test rad2deg(fa.λ_Sa) ≈ 102.168400 atol=1e-7 rtol=1e-7
        # Testing Mean Longitude of Uranus
        @test rad2deg(fa.λ_Ur) ≈ 332.317825 atol=1e-7 rtol=1e-7
        # Testing Mean Longitude of Neptune (from SOFA routines)
        @test fa.λ_Ne ≈ 5.474423134426369 atol=1e-20 rtol=1e-20
        # Testing General Accumulated Precession 
        @test rad2deg(fa.pₐ) ≈ 0.059545 atol=1e-6 rtol=1e-7
    end 

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

