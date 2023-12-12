
@testset "legacy" verbose = true begin

    # Test IDs 
    @test Orient.AXESID_ECLIPB1950 == 18
    @test Orient.AXESID_FK4 == 3
    @test Orient.AXESID_GALACTIC == 13
    @test Orient.AXESID_B1950 == 2

    # Function to compute the angle between 2 vectors in arcseconds
    v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600

    v = rand(BigFloat, 3)
    v /= norm(v)

    R = Orient.DCM_MEME2000_TO_B1950
    Rₑ = pxform("J2000", "B1950", 0.0)
    @test v2as(R*v, Rₑ*v) ≤ 1e-10

    R = Orient.DCM_B1950_TO_FK4
    Rₑ = pxform("B1950", "FK4", 0.0)
    @test v2as(R*v, Rₑ*v) ≤ 1e-10

    R = Orient.DCM_FK4_TO_GALACTIC
    Rₑ = pxform("FK4", "GALACTIC", 0.0)
    @test v2as(R*v, Rₑ*v) ≤ 1e-10

    R = Orient.DCM_B1950_TO_ECLIPB1950
    Rₑ = pxform("B1950", "ECLIPB1950", 0.0)
    @test v2as(R*v, Rₑ*v) ≤ 1e-10

end;
