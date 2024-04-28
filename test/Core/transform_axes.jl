using Test
using LinearAlgebra
using FrameTransformations
using ReferenceFrameRotations
using Tempo

g = FrameSystem{2, Float64, BarycentricDynamicalTime}()

@test rotation3(g, 0, 0, 1.0) == Rotation{1, Float64}(1.0I)
@test rotation6(g, 0, 0, 1.0) == Rotation{2, Float64}(1.0I)
# order is to low
@test_throws Exception rotation9(g, 0, 0, 1.0) 
@test_throws Exception rotation12(g, 0, 0, 1.0) 

g = FrameSystem{4, Float64, BarycentricDynamicalTime}()

@test rotation9(g, 0, 0, 1.0) == Rotation{3, Float64}(1.0I)
@test rotation12(g, 0, 0, 1.0) == Rotation{4, Float64}(1.0I)

# No node in the frame system
@test_throws Exception rotation3(g, 0, 1, 1.0)
@test_throws Exception rotation6(g, 0, 1, 1.0)
@test_throws Exception rotation9(g, 0, 1, 1.0)
@test_throws Exception rotation12(g, 0, 1, 1.0)

# ---
# Root axes

faxs0 = FrameTransformations.FrameAxesFunctions{4, Float64}()
fax0 = FrameTransformations.FrameAxesNode{4, Float64}(:Ax0, 0, 1, 1, faxs0)

@test fax0.name == :Ax0
@test length(faxs0.fun) == 4 
@test fax0.class == 0 
@test fax0.id == 1
@test fax0.parentid == 1
@test fax0.f === faxs0

# ---
# Direct rotation
f(t) = Rotation{4}(angle_to_dcm(π/4, :Z))
faxs1 = FrameTransformations.FrameAxesFunctions{4, Float64}(f)

@test faxs1[1](0.0) === f(0.0)

fax1 = FrameTransformations.FrameAxesNode{4, Float64}(:Ax1, 0, 2, 1, faxs1)

@test fax1.name == :Ax1
@test fax1.class == 0 
@test fax1.id == 2
@test fax1.parentid == 1
@test fax1.f === faxs1

# ---
# Inv rotation
finv(t) = inv(f(t))
faxs2 = FrameTransformations.FrameAxesFunctions{4, Float64}(finv)
fax2 = FrameTransformations.FrameAxesNode{4, Float64}(:Ax2, 0, 3, 2, faxs2)

for ax in [fax0, fax1, fax2]
    add_axes!(g, ax) 
    FrameTransformations.add_edge!(g.axes, ax.parentid, ax.id)
end

@test rotation3(g, 1, 3, 0.0)[1] ≈ DCM(I)
@test FrameTransformations.order(rotation3(g, 1, 1, 0.0)) == 1
@test FrameTransformations.order(rotation6(g, 1, 1, 0.0)) == 2
@test FrameTransformations.order(rotation9(g, 1, 1, 0.0)) == 3
@test FrameTransformations.order(rotation12(g, 1, 1, 0.0)) == 4


# Dummy constructor 
funs = FrameTransformations.FrameAxesFunctions{4, Int}()
@test typeof(funs) == FrameTransformations.FrameAxesFunctions{4, Int, 12}
@test length(funs.fun) == 4

for i in 1:3
    @test funs.fun[i](0) == Rotation{4}(1I)
end

fcn_rot(t) = Rotation{2}(angle_to_dcm(t, :Y))

# Default constructor
funs = FrameTransformations.FrameAxesFunctions{Float64}(fcn_rot, fcn_rot)
@test typeof(funs) == FrameTransformations.FrameAxesFunctions{2, Float64, 6}
@test length(funs.fun) == 2

@test funs[2](π / 3)[1] ≈ fcn_rot2(π / 3)[1] atol = 1e-12
@test funs[2](π / 3)[2] ≈ fcn_rot2(π / 3)[2] atol = 1e-12
