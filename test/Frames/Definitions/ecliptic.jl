
# create dummy frame system 
FRAMES = FrameSystem{3, Float64}()

# insert ICRF
@axes ICRF 1 InternationalCelestialReferenceFrame
@axes MEME2000 2 
@axes ECLIPJ2000 3

add_axes_inertial!(FRAMES, ICRF)
add_axes_meme2000!(FRAMES, MEME2000, ICRF)
add_axes_eclipj2000!(FRAMES, ECLIPJ2000, MEME2000)

@testset "Ecliptic Equinox at J2000" begin
    ep = rand(0.0:1e7)

    R = sxform("J2000", "ECLIPJ2000", ep)
    R_ = rotation6(FRAMES, MEME2000, ECLIPJ2000, ep)

    @test isapprox(R[1:3, 1:3],  R_[1], atol=1e-8)

end