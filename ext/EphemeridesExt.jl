module EphemeridesExt

import FrameTransformations: add_point_ephemeris!, add_axes_ephemeris!

using FrameTransformations: FrameSystem,
    FramePointFunctions, Translation, add_point!,
    FrameAxesFunctions, Rotation, add_axes!,
    check_point_ephemeris, check_axes_ephemeris,
    angles_to_rot3, angles_to_rot6, angles_to_rot9, angles_to_rot12

using Ephemerides: EphemerisProvider,
    ephem_vector3, ephem_vector6, ephem_vector9, ephem_vector12,
    ephem_rotation3, ephem_rotation6, ephem_rotation9, ephem_rotation12


"""
    add_point_ephemeris!(fr::FrameSystem{O, T}, eph::EphemerisProvider, 
        name::Symbol, id::Int) where {O, T}
    
Add a point from `Ephemerides.jl` provider.
"""
function add_point_ephemeris!(
    fr::FrameSystem{O,T}, eph::EphemerisProvider, name::Symbol, id::Int
) where {O,T}

    pid, axid = check_point_ephemeris(fr, eph, id)

    funs = FramePointFunctions{O,T}(
        t -> Translation{O}(ephem_vector3(eph, pid, id, t)),
        t -> Translation{O}(ephem_vector6(eph, pid, id, t)),
        t -> Translation{O}(ephem_vector9(eph, pid, id, t)),
        t -> Translation{O}(ephem_vector12(eph, pid, id, t))
    )
    return add_point!(fr, name, id, axid, funs, pid)
end

"""
    add_axes_ephemeris!(fr::FrameSystem{O, T}, eph::EphemerisProvider, 
        name::Symbol, id::Int) where {O, T}
    
Add an axes from `Ephemerides.jl` provider.
"""
function add_axes_ephemeris!(
    fr::FrameSystem{O,T}, eph::EphemerisProvider, name::Symbol, id::Int, rot_seq::Symbol
) where {O,T}

    # Check and retrieve the parent ID for the given axes
    pid = check_axes_ephemeris(fr, eph, id)

    if rot_seq in (:ZYX, :XYX, :XYZ, :XZX, :XZY, :YXY, :YXZ, :YZX, :YZY, :ZXY, :ZXZ, :ZYZ)
        funs = FrameAxesFunctions{O,T}(
            t -> Rotation{O}(angles_to_rot3(ephem_rotation3(eph, pid, id, t), rot_seq)),
            t -> Rotation{O}(angles_to_rot6(ephem_rotation6(eph, pid, id, t), rot_seq)),
            t -> Rotation{O}(angles_to_rot9(ephem_rotation9(eph, pid, id, t), rot_seq)),
            t -> Rotation{O}(angles_to_rot12(ephem_rotation12(eph, pid, id, t), rot_seq))
        )
    else
        throw(ArgumentError("The rotation sequence :$rot_seq is not valid."))
    end
    return add_axes!(fr, name, id, funs, pid)
end

end