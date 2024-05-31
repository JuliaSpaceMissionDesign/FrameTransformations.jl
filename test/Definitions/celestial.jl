using Test
using LinearAlgebra
using FrameTransformations
using ReferenceFrameRotations

frames = FrameSystem{3,Float64}()

# Add ICRF axes
@test_nowarn add_axes_root!(frames, :ICRF, 1)

# test that ICRF can only be added as a root axes
@test_throws ArgumentError add_axes_icrf!(frames)

# Test GCRF can be added as child of ICRF
@test_nowarn add_axes_gcrf!(frames)

frames = FrameSystem{3, Float64}() 
add_axes_icrf!(frames)

node = axes_graph(frames).nodes[1]
@test node.id == AXESID_ICRF
@test node.name == :ICRF

add_axes_gcrf!(frames)
node = axes_graph(frames).nodes[2]
@test node.id == AXESID_GCRF
@test node.name == :GCRF 

# Test they have an identity rotation
R = rotation6(frames, 1, 23, 0.0)
@test R[1] ≈ DCM(1.0I) atol=1e-12 rtol=1e-12
@test R[2] ≈ DCM(0.0I) atol=1e-12 rtol=1e-12

# Test GCRF as root axes 
frames = FrameSystem{3, Float64}() 
add_axes_gcrf!(frames)

node = axes_graph(frames).nodes[1]
@test node.id == AXESID_GCRF
@test node.name == :GCRF 
