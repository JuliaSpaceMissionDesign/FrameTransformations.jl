# # [Use Case: CR3BP](@id example_01_cr3bp)
# _This example was generated on DATEOFTODAY._

# The power of the [`FrameSystem`](@ref) is its capability to handle axes 
# transformations and point translations of both high-accuracy and simplified 
# models. The use-case here presented includes the case of the Circular-Restricted 
# Three-Body Problem (CR3BP) rotating frame transformation handling.

# In particular, when dealing with the [CR3BP](https://orbital-mechanics.space/the-n-body-problem/circular-restricted-three-body-problem.html), 
# mission analysis are used to exploit non-dimensional, rotating coordinates 
# to express the equations of motion and perform the computations. 

# In this tutorial, we create a [`FrameSystem`](@ref) to handle transformations 
# within the Earth-Moon CR3BP, which is characterized by a mass ratio of 
# approximately `μ = 0.012`. We start off by creating a frame system without any 
# ephemeris provider, since we are using a simplified model.

using FrameTransformations

CR3BP = FrameSystem{2, Float64}()

# As always, the first step requires the definition of the root axes and points. 
# In this case, we use the a generic set of _inertial axes_ and the Earth-Moon 
# Barycenter (EMB).

@axes InertialAx 1 InertialFrame 

add_axes_inertial!(CR3BP, InertialAx)

@point EMBc 1 EarthMoonBarycenterCr3bp

add_point_root!(CR3BP, EMBc, InertialAx)

# We now proceed to add our synodic axes: in the CR3BP these are uniformly 
# rotating with respect to the `InertialAx` about the Z-axis. Therefore, we 
# leverage the [rotating axes](@ref rot_axes) type: 

using ReferenceFrameRotations

f(t) = angle_to_dcm(t, :Z) 

@axes SynodicAx 2 SynodicFrame 

add_axes_rotating!(CR3BP, SynodicAx, InertialAx, f)

# Note that there is no need to specify the rotation derivatives, as they'll 
# be computed by  automatic differentiation via the 
# [ForwardDiff](https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl) package. 
# For performace-critical transformations, however, it is reccomended to manually 
# define these derivatives.

# Now, let's assume we have our spacecraft. Most likely, its states will be 
# expressed in the synodic frame. In this case, we leverage 
# [updatable points](@ref updatable_points), since we desire to manually update 
# its state at each time.

@point Spacecraft -1_900_000

add_point_updatable!(CR3BP, Spacecraft, EMBc, SynodicAx)

# We know assume that at `t = 0.8` our spacecraft is at L4, therefore we update 
# its state accordingly:

μ = 0.012
xL4 = [1/2-μ, sqrt(3)/2, 0.0, 0.0, 0.0, 0.0]
t = 0.8

update_point!(CR3BP, Spacecraft, xL4, t)

# Finally we can retrieve the spacecraft state in both the synodic as well as the 
# inertial axes: 

vector6(CR3BP, EMBc, Spacecraft, SynodicAx, t)

#- 
vector6(CR3BP, EMBc, Spacecraft, InertialAx, t)