
@axes ICRF 1 InternationalCelestialReferenceFrame
@axes GCRF 23
@axes TestAxes 2

@testset "Celestial" verbose=true begin 

    v2as = (x, y) -> acosd(max(-1, min(1, dot(x / norm(x), y / norm(y))))) * 3600

    frames = FrameSystem{3,Float64}()
    add_axes_inertial!(frames, TestAxes)

    # test that ICRF can only be added as a root axes
    @test_throws ArgumentError add_axes_icrf!(frames)

    # test that GCRF can only be added as root or parent of ICRF 
    @test_throws ArgumentError add_axes_gcrf!(frames)
    
    frames = FrameSystem{3, Float64}() 
    add_axes_icrf!(frames)

    node = frames_axes(frames).nodes[1]
    @test node.id == AXESID_ICRF
    @test node.name == :ICRF

    add_axes_gcrf!(frames)
    node = frames_axes(frames).nodes[2]
    @test node.id == AXESID_GCRF
    @test node.name == :GCRF 

    # Test they have an identity rotation
    R = rotation6(frames, ICRF, GCRF, 0.0)
    @test R[1] ≈ DCM(1.0I) atol=1e-12 rtol=1e-12
    @test R[2] ≈ DCM(0.0I) atol=1e-12 rtol=1e-12

    # Test GCRF as root axes 
    frames = FrameSystem{3, Float64}() 
    add_axes_gcrf!(frames)

    node = frames_axes(frames).nodes[1]
    @test node.id == AXESID_GCRF
    @test node.name == :GCRF 

end;