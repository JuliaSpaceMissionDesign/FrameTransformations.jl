
# create dummy frame system 
FRAMES = FrameSystem{3, Float64}()

# insert ICRF
@axes ICRF 1 InternationalCelestialReferenceFrame
@axes MEME2000 22 
@axes ECLIPJ2000 17

add_axes_inertial!(FRAMES, ICRF)
add_axes_meme2000!(FRAMES, MEME2000, ICRF)
add_axes_eclipj2000!(FRAMES, ECLIPJ2000, MEME2000)

@testset "Ecliptic Equinox at J2000" begin
    ep = rand(0.0:1e7)

    v2as = (x, y) -> acosd(max(-1, min(1, dot(x/norm(x), y/norm(y)))))*3600
    
    R = sxform("J2000", "ECLIPJ2000", ep)
    R_ = rotation6(FRAMES, MEME2000, ECLIPJ2000, ep)
    
    v = rand(BigFloat, 3); 
    v /= norm(v)

    @test v2as(R[1:3, 1:3]*v, R_[1]*v) ≈ 0.0 atol=1e-6
    @test R_[2] ≈ zeros(3,3) atol=1e-14

end;