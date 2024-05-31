module FrameTransformations 

using LinearAlgebra
import LinearAlgebra: UniformScaling, matprod

using ForwardDiff
using StaticArrays
using ReferenceFrameRotations
using FunctionWrappers: FunctionWrapper
using FunctionWrappersWrappers: FunctionWrappersWrapper

using JSMDUtils
using JSMDUtils.Autodiff
using JSMDUtils.Math: D¹, D², D³, 
                      arcsec2rad,
                      unitvec, δunitvec, δ²unitvec, δ³unitvec, 
                      cross3, cross6, cross9, cross12, 
                      _3angles_to_δdcm, _3angles_to_δ²dcm, _3angles_to_δ³dcm

using JSMDInterfaces.Graph: AbstractJSMDGraphNode, add_edge!, add_vertex!, get_path, has_vertex

using JSMDInterfaces.Ephemeris: AbstractEphemerisProvider, 
                                ephem_position_records, ephem_available_points,
                                ephem_orient_records, ephem_available_axes
using JSMDInterfaces.Interface: @interface

using IERSConventions: iers_bias, iers_obliquity, 
                       iers_rot3_gcrf_to_itrf, iers_rot6_gcrf_to_itrf, 
                       iers_rot9_gcrf_to_itrf, iers_rot12_gcrf_to_itrf,
                       iers_rot3_gcrf_to_mod, iers_rot3_gcrf_to_tod, iers_rot3_gcrf_to_gtod, 
                       iers_rot3_gcrf_to_pef, iers_rot3_gcrf_to_cirf, iers_rot3_gcrf_to_tirf,
                       IERSModel, iers2010a, iers2010b, iers1996

using SMDGraphs: MappedNodeGraph, SimpleGraph, MappedGraph,
                 get_mappedid, get_mappednode,  get_node, get_path

import SMDGraphs: get_node_id

using Tempo: AbstractTimeScale, Epoch, j2000s, BarycentricDynamicalTime

# Core Routines
# ==============================================

# Low-level types and aliases
export Rotation

include("Core/rotation.jl")
include("Core/ad.jl")

# Frame system 
export FrameSystem, 
       order, timescale, points_graph, axes_graph, directions_map, 
       has_axes, has_point, has_direction,
       point_id, axes_id,
       axes, points, directions,
       rotation3, rotation6, rotation9, rotation12, 
       vector3, vector6, vector9, vector12,
       direction3, direction6, direction9, direction12

include("Core/nodes.jl")
include("Core/graph.jl")
include("Core/transform.jl")

# Helper functions 
export add_axes!, add_point!, add_direction!,
       add_axes_root!, add_axes_inertial!, add_axes_rotating!, add_axes_fixedoffset!, add_axes_alias!,
       add_point_root!, add_point_dynamical!, add_point_fixedoffset!, add_point_alias!,
       add_direction_fixed!, add_direction_position!, add_direction_velocity!, 
       add_direction_orthogonal!, add_direction_normalize!

include("Core/axes.jl")
include("Core/points.jl")
include("Core/directions.jl")

# Standard axes/points definitions
# ==============================================

export AXESID_ICRF, AXESID_GCRF, 
       AXESID_ECL2000, AXESID_EME2000, 
       AXESID_MOONME_DE421, AXESID_MOONPA_DE421, AXESID_MOONPA_DE440

include("Definitions/index.jl")

export add_point_ephemeris!, add_axes_ephemeris!

include("Definitions/ephemeris.jl")

export add_axes_icrf!, add_axes_gcrf!, add_axes_eme2000!, add_axes_ecl2000!

include("Definitions/celestial.jl")
include("Definitions/ecliptic.jl")

export add_axes_itrf!, add_axes_cirf!, add_axes_tirf!,
       add_axes_mod!, add_axes_tod!, add_axes_gtod!, add_axes_pef!

include("Definitions/terrestrial.jl")

export add_axes_pa440!, add_axes_pa421!, add_axes_me421!

include("Definitions/lunar.jl")

export add_point_surface!, add_axes_topocentric!

include("Definitions/topocentric.jl")

export  add_axes_twovectors!

include("Definitions/twovector.jl")

end