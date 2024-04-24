module EphemeridesExt 

import FrameTransformations: add_point_ephemeris!, add_axes_ephemeris!

using FrameTransformations: FrameSystem, 
                            FramePointFunctions, SVectorNT, add_point!,
                            FrameAxesFunctions, Rotation, add_axes!, 
                            triplet_to_rot3, triplet_to_rot6, triplet_to_rot9, triplet_to_rot12,
                            check_point_ephemeris, check_axes_ephemeris, 
                            POINT_CLASSID_DYNAMIC

using Ephemerides: EphemerisProvider, 
                   ephem_vector3, ephem_vector6, ephem_vector9, ephem_vector12,
                   ephem_rotation3, ephem_rotation6, ephem_rotation9, ephem_rotation12 

function add_point_ephemeris!(
    frames::FrameSystem{O, N}, eph::EphemerisProvider, name::Symbol, id::Int
) where {O, N}

    parentid, axesid = check_point_ephemeris(frames, eph, id)

    funs = FramePointFunctions{O, N}(
        t -> SVectorNT{3O, N}( ephem_vector3(eph, parentid, id, t) ),
        t -> SVectorNT{3O, N}( ephem_vector6(eph, parentid, id, t) ),
        t -> SVectorNT{3O, N}( ephem_vector9(eph, parentid, id, t) ),
        t -> SVectorNT{3O, N}( ephem_vector12(eph, parentid, id, t) )
    )

    return add_point!(frames, name, id, axesid, POINT_CLASSID_DYNAMIC, funs, parentid)    
end

function add_axes_ephemeris!(
    frames::FrameSystem{O,T}, eph::EphemerisProvider, name::Symbol, id::Int, class::Int, rot_seq::Symbol
) where {O,T}

    # Check and retrieve the parent ID for the given axes
    parentid = check_axes_ephemeris(frames, id)

    if rot_seq in (:ZYX, :XYX, :XYZ, :XZX, :XZY, :YXY, :YXZ, :YZX, :YZY, :ZXY, :ZXZ, :ZYZ)
        funs = FrameAxesFunctions{T,O}(
            t -> Rotation{O}(triplet_to_rot3(ephem_rotation3(eph, parentid, id, t), rot_seq)),
            t -> Rotation{O}(triplet_to_rot6(ephem_rotation6(eph, parentid, id, t), rot_seq)),
            t -> Rotation{O}(triplet_to_rot9(ephem_rotation9(eph, parentid, id, t), rot_seq)),
            t -> Rotation{O}(triplet_to_rot12(ephem_rotation12(eph, parentid, id, t), rot_seq))
        )
    else
        throw(ArgumentError("The rotation sequence :$rot_seq is not valid."))
    end

    return add_axes!(frames, name, id, class, funs, parentid)
end

end