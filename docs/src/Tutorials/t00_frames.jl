# # [Frame System Overview](@id tutorial_00_frames)
# _This example was generated on DATEOFTODAY._

# The core object of `FrameTransformations` is the [`FrameSystem`](@ref), which provides 
# the capability to compute relative position, orientation and their time derivatives up to 
# order 3 (jerk), between standard and user-defined point and axes. It works by creating two 
# separate graphs that silently store and manage all the parent-child relationships between 
# the user-registered axes and points, in the form of `FramePointNode` and `FrameAxesNode`. 

# These two objects define two precise entities: 
# - **Axes**: defines an orientation in space. These are related each other by means of a 
#   [`Rotation`](@ref) transformation which relate one axes to a parent axes in 
#   a certain time interval.
# - **Points**: defines a location in space. These are related each other by 
#   means of a `Translation` transformation which relate one point to a parent point in a 
#   particular axes in a certain time interval.

# Additionally, it is possible to create `Direction`s, as vector valued functions that could
# be used to define custom frames.

#-
#md # !!! note 
#md #     A single [`FrameSystem`](@ref) instance simultaneously handles both the axes and 
#md #     point graphs, regardless of what the user registers in it. For instance, if no 
#md #     points are added, the point graph will remain empty. The same applies for directions.
#-

# Additionally, any node can have several childs, each with different transformations with 
# respect to the parent node. However, they shall be **registered** within the 
# [`FrameSystem`](@ref) before being used in a transformation or as parents of other nodes.


# ## Basic Constructors 
# The creation of a generic [`FrameSystem`](@ref) requires the definition of the maximum 
# desired transformation order and of its `DataType`, which in most applications is a `Float64`. 
# The transformation order is always one greater than the maximum desired time derivative.
# For instance, if the user only desires to compute position and velocity components (i.e., 
# order 1 time-derivative), the transformation order to be used is 2. Thus, the maximum 
# allowed transformation order is 4. 

# In this example, we highlight the most basic way to initialise a [`FrameSystem`](@ref):

using FrameTransformations
using Tempo 

F = FrameSystem{2, Float64}()

# From this example, you can see that within the frame system there are both point and axes 
# graphs. However, at the moment they are completely empty since the graph was just created.

# Each [`FrameSystem`](@ref) object is assigned a reference timescale that is used to perform 
# computations with epochs and to parse ephemeris files. The default timescale is the 
# `BarycentricDynamicalTime`, however, the user is free to select the most suited timescale 
# for his applications. In this example, we set the `InternationalAtomicTime` as the reference scale.

F = FrameSystem{2, Float64, InternationalAtomicTime}()

# ## Graph Inspection

# Once a [`FrameSystem`](@ref) is constructed (and populated) there are many routines devoted 
# to inspect its content. As already said, there are three main *objects* that are contained 
# in the `FrameSystem`: **points**, **axes** and **directions**. For each of them series of 
# utility functions are made available in order to check for the presence of a registered point:

has_point(F, 1)

# a registered axes:

has_axes(F, 1)

# or a registered direction:

has_direction(F, :Root)

# Additionally, the possibility to get a dictionary containing all name-id relationships is 
# made available for axes, via the [`axes`](@ref) method:

axes(F)

# and points, via the [`points`](@ref) method:

points(F)

# Finally, the `FrameSystem` order and timescale might be retrieved via the associated methods:

order(F)

#- 

FrameTransformations.timescale(F)

# Refer to the [API](@ref frames_api) for additional details.

# ## Basic Usage

#md # !!! note 
#md #     Work in progress

# ## Ephemerides Support

# In certain scenarios, the transformations require usage of binary ephemeris kernels, e.g., 
# the JPL's DE440 files. To support this applications, this package has an interface relying 
# on [JSMDInterfaces.jl](https://github.com/JuliaSpaceMissionDesign/JSMDInterfaces.jl) 
# `AbstractEphemerisProvider`s. Currently, this package is shipped with extension for the 
# following two ephemeris readers:
# * [Ephemerides.jl](https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl)
# * [CalcephEphemeris.jl](https://github.com/JuliaSpaceMissionDesign/CalcephEphemeris.jl)

# Once the desired ephemeris provider is created, it can be used to register points or axes. 
# In this example we begin loading an old DE421 kernels to pass to the ephemeris reader.

using Ephemerides, Downloads

url = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp";
eph = EphemerisProvider(Downloads.download(url));

F = FrameSystem{2, Float64}()

# Before registering any node, a set of root axes and a root node shall be anyway registered.

add_axes_icrf!(F)
add_point_root!(F, :SSB, 0, 1)

# Points from the `EphemerisProvider` can be now registered. 

add_point_ephemeris!(F, eph, :Sun, 10)
add_point_ephemeris!(F, eph, :EMB, 3)

# Here the parent point will be inferred from the ephemeris.