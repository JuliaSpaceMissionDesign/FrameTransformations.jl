# # [Points Creation and Translations](@id tutorial_02_points)
# _This example was generated on DATEOFTODAY._

# Similarly to [axes](@ref tutorial_01_axes), `FrameTransformations` also provides the 
# capability to define custom and standard reference points (e.g., the Solar System 
# Barycenter) and arbitrarily connect them through the [`FrameSystem`](@ref). In turn, this 
# allows the computation of the relative position and its derivatives (up to order 3) between 
# any two registered points and express it in any known set of axes. 

# At the time being, the following types of points are supported:
# - **Root point**: it is the root of the point graph.
# - **Fixed points**: are those whose positions have a constant offset with respect their 
#   parent point in a given set of axes.
# - **Dynamical points**: the position of these points depends only on time and is computed 
#   through custom user-defined functions.
# - **Ephemeris points**: are those whose state-vector is retrieved from binary SPK kernels 
#   (e.g., DE440) that are loaded within the [`FrameSystem`](@ref).

#md # !!! note
#md #     This package provides a dedicated function to register each type of supported points.

# ## Graph Initialisation

# In this section we will display how to create a frame system to compute generic points 
# transformation. Differently from the axes graph, each register point is also associated 
# to a set of axes. Hence, this tutorial assumes the reader to already be familiar with the 
# different types of axes and their definitions. 

# We then can go ahead and initialise the graph.

using StaticArrays
using FrameTransformations

F = FrameSystem{2,Float64}()

# ## Root Point

# To initialise the point graph, we first need to define a root point. This, in turn, must 
# be associated to an arbitrary set of axes. Therefore, we begin by definining a generic 
# `SatFrame`, here considered as inertial, and then register a root point, called 
# `SC` in our graph. 

# A root point can be registered using the [`add_point!`](@ref) function: 

add_axes!(F, :SatFrame, -1)

add_point!(F, :SC, -10000, :SatFrame)

#md # !!! tip 
#md #     For standard applications, it is good practice that the points's IDs are as in 
#md #     agreement with NAIF's numbering system. This becomes mandatory to properly read 
#md #     JPL's SPK kernels.

#md # !!! note 
#md #     The frame system uses an integer system based on the user-defined IDs to compute 
#md #     the transformations between axes and points. The name and acronym of the point are 
#md #     only used as aliases to provide a user-friendly interface to the transformations 
#md #     and do not have any other meaning. 

# ## Fixed Points

# Fixed points have a constant relative position vector with respect to their parent points 
# in a given set of axes. Similarly to fixed-offset axes, these points are fixed w.r.t. their 
# parents but might be moving with respect to others. 

# In this example, we use the [`add_point_fixedoffset!`](@ref) function to register the location 
# of an antenna and two solar panels, which are generally fixed in the satellite body-fixed 
# frame. To do so, we define a position offset in the form of a 3-elements vector with respect 
# to the `SC`.

sa_offset_left = [1.0, 0.0, 0.0]
sa_offset_right = [-1.0, 0.0, 0.0]
an_offset = [0.0, 0.0, -1.0]

add_point_fixedoffset!(F, :SolArrLeft, -10101, :SC, :SatFrame, sa_offset_left)
add_point_fixedoffset!(F, :SolArrRight, -10102, :SC, :SatFrame, sa_offset_right)
add_point_fixedoffset!(F, :Antenna, -10001, :SC, :SatFrame, an_offset)

# As a result the graph is now populated with the new points and we can finally compute 
# their relative positions and velocities with the proper transformation functions: 

#-
vector3(F, :SolArrLeft, :SC, :SatFrame, 123.0)

#-
vector6(F, :Antenna, :SolArrRight, :SatFrame, 456.0)

# As expected, since these points are fixed, the relative velocity vector is null.

# ## Dynamical Points

# Dynamical points are generic time-dependent points whose position vector (and optionally 
# its derivatives) are only function of time. However, differently from ephemeris points, 
# their position is computed through user-defined functions.

fun(t) = SA[cos(t), sin(t), 0.0]

add_point_dynamical!(F, :TimedAppendage, -10003, :SolArrLeft, :SatFrame, fun)

#-
vector3(F, :TimedAppendage, :SC, :SatFrame, π / 3)

#md # !!! note 
#md #     To avoid allocations, `fun` should return a static array.

# Similarly to rotating-axes, if the user only provides the function to compute the relative 
# position, the remaining derivatives are automatically retrievied via automatic 
# differentiation of `fun`. On the other hand, if those functions are specified, they must 
# return a single vector that stacks all the components. For instance, for the first order 
# derivative of `fun`, the function should return a 6-elements vector containing the
# relative position and velocity. For example: 

fun(t) = SA[cos(t), sin(t), 0]
dfun(t) = SA[cos(t), sin(t), 0, -sin(t), cos(t), 0]

add_point_dynamical!(F, :TimedAppendage2, -10004, :SolArrLeft, :SatFrame, fun, dfun)

#- 
vector6(F, :TimedAppendage2, :SC, :SatFrame, π / 3)

# We can again see that the results are in agreement with the previous example. 
# For more details, consult the [`add_point_dynamical!`](@ref) documentation.

# ## Ephemeris Points

# Refer to the [frames tutorial](@ref tutorial_00_frames)'s *Ephemeris Support* section.