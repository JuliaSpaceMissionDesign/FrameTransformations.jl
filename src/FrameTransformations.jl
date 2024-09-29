module FrameTransformations

using LinearAlgebra
using StaticArrays
using ReferenceFrameRotations
using FunctionWrappers: FunctionWrapper
using FunctionWrappersWrappers: FunctionWrappersWrapper

using JSMDUtils.Math: D¹, D², D³, arcsec2rad,
       unitvec, δunitvec, δ²unitvec, δ³unitvec,
       cross3, cross6, cross9, cross12,
       angle_to_δdcm, angle_to_δdcm, angle_to_δ²dcm, angle_to_δ³dcm

using JSMDInterfaces.Graph: AbstractJSMDGraphNode,
       add_edge!, add_vertex!, get_path, has_vertex

using SMDGraphs: MappedNodeGraph, SimpleGraph, MappedGraph,
       get_mappedid, get_mappednode, get_node, get_path

import SMDGraphs: get_node_id

using Tempo: AbstractTimeScale, Epoch, j2000s, BarycentricDynamicalTime, ftype

using JSMDInterfaces.Ephemeris: AbstractEphemerisProvider,
       ephem_position_records, ephem_available_points,
       ephem_orient_records, ephem_available_axes

using JSMDInterfaces.Interface: @interface

using JSMDInterfaces.Bodies: body_rotational_elements, ∂body_rotational_elements,
       ∂²body_rotational_elements, ∂³body_rotational_elements

using IERSConventions: iers_bias, iers_obliquity,
       iers_rot3_gcrf_to_itrf, iers_rot6_gcrf_to_itrf,
       iers_rot9_gcrf_to_itrf, iers_rot12_gcrf_to_itrf,
       iers_rot3_gcrf_to_mod, iers_rot3_gcrf_to_tod, iers_rot3_gcrf_to_gtod,
       iers_rot3_gcrf_to_pef, iers_rot3_gcrf_to_cirf, iers_rot3_gcrf_to_tirf,
       IERSModel, iers2010a, iers2010b, iers1996

using ForwardDiff
using JSMDUtils.Autodiff: JSMDDiffTag, derivative

# ==========================================================================================
# Core
# ==========================================================================================

# Low-level types and aliases
export Translation, Rotation

include("Core/translation.jl")
include("Core/rotation.jl")
include("Core/ad.jl")

# Frame system 
export FrameSystem,
       order, timescale, points_graph, axes_graph, points_alias, axes_alias, directions,
       has_axes, has_point, has_direction,
       point_id, axes_id

include("Core/nodes.jl")
include("Core/graph.jl")

# Helper functions 
export add_axes!, add_axes_projected!, add_axes_rotating!, add_axes_fixedoffset!,
       add_point!, add_point_dynamical!, add_point_fixedoffset!,
       add_direction!

include("Core/axes.jl")
include("Core/points.jl")
include("Core/directions.jl")

# Transformations 
export rotation3, rotation6, rotation9, rotation12,
       vector3, vector6, vector9, vector12,
       direction3, direction6, direction9, direction12

include("Core/transform.jl")

# ==========================================================================================
# Definitions
# ==========================================================================================

export AXESID_ICRF, AXESID_GCRF,
       AXESID_ECL2000, AXESID_EME2000,
       AXESID_MOONME_DE421, AXESID_MOONPA_DE421, AXESID_MOONPA_DE440

include("Definitions/index.jl")

export add_axes_icrf!, add_axes_gcrf!, add_axes_eme2000!, add_axes_ecl2000!

include("Definitions/celestial.jl")
include("Definitions/ecliptic.jl")

export add_point_ephemeris!

include("Definitions/ephemeris.jl")

export add_axes_frozen!

include("Definitions/frozen.jl")

export add_axes_itrf!, add_axes_cirf!, add_axes_tirf!,
       add_axes_mod!, add_axes_tod!, add_axes_gtod!, add_axes_pef!

include("Definitions/terrestrial.jl")

export add_axes_bci2000!, add_axes_bcrtod!

include("Definitions/planetary.jl")

export add_axes_pa440!, add_axes_pa421!, add_axes_me421!

include("Definitions/lunar.jl")

export add_axes_topocentric!, add_point_surface!

include("Definitions/topocentric.jl")

export add_direction_position!, add_direction_velocity!, add_direction_orthogonal!,
       add_direction_fixed!

include("Definitions/directions.jl")

export add_axes_twodir!

include("Definitions/axesfromdir.jl")

export add_axes_fixed_quaternion!, add_axes_fixed_angles!, add_axes_fixed_angleaxis!

include("Definitions/attitude.jl")

end