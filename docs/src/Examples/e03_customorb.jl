# # [Use Case: Custom spacecraft orbit](@id example_03_orb)
# _This example was generated on DATEOFTODAY._

# Once the general structure of the [`FrameSystem`](@ref) is understood, we can pass to 
# a use case in which we want to build and exploit our frame system to perform computations 
# inserting a custom orbit. This becomes especially crucial in complex cases like Trajectory 
# Optimization and Navigation Analysis.

# In such scenarios, the trajectory is under design, and the trajectory information might 
# not be completely available. Moreover, in these cases derivatives could be required for 
# various quantities such as time, states, and parameters.

# In this context, we shall remember that `FrameTransformations` is able to perform operations,
# including AD, on the frames **independent variable**, e.g. only time.

# A proper orbit representation is essential to avoid perturbation confusion and ensure 
# proper custom orbit handling. For this purpose, two point types can seems suitable: 
# `updatable` and `dynamical` points. 
# In this case, however, `updatable` points are not well-suited as they are essentially 
# constants for the AD system. Then, `dynamical` points can effectively handle this scenario. 

# ### Frame system setup

# First of all, a new [`FrameSystem`](@ref) shall be created. 

using FrameTransformations 

G = FrameSystem{4, Float64}()

add_axes_icrf!(G)
add_point_root!(G, :Earth, 399, 1)

#md # !!! note 
#md #     As the orbit considered in this example is custom, there is no requirement 
#md #     to load ephemeris. 

# ### Custom orbit model

# A custom orbit model is then required. In this case two dummy orbits are 
# created and stored in a `CubicSpline` object.

using JSMDInterfaces.Math: interpolate
using JSMDUtils.Math: InterpCubicSplines

c = 2π/86400
de = 0:3600:2*86400;
fv(t) = [cos(c*t), sin(c*t), 0, -c*sin(c*t), c*cos(c*t), 0];
fv2(t) = [sin(c*t), cos(c*t), 0];
const ORB1 = InterpCubicSplines(de, hcat([fv(e) for e in de]...));
const ORB2 = InterpCubicSplines(de, hcat([fv2(e) for e in de]...));

# ### Custom orbit handling 

# Given the custom orbit model, [`add_point_dynamical!`](@ref) could be used 
# to link it to the frame system.

# Insert dynamical points
add_point_dynamical!(G, :SC1, -1, 399, 1, t->interpolate(ORB1, t)[1:3], t->interpolate(ORB1, t))
add_point_dynamical!(G, :SC2, -2, 399, 1, t->interpolate(ORB2, t))

# Once in the frame system, AD could be exploited to retrieve states+gradients w.r.t. 
# time in any registered frame:
 
vector12(G, 399, -1, 1, 12345.0)

# It is also possible to use the computational graph to compute quantities between 
# completely custom states representations:

vector12(G, -1, -2, 1, 12345.0)