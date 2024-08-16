module EphemeridesExt 

import FrameTransformations: add_point_ephemeris!

using FrameTransformations: FrameSystem, 
                            FramePointFunctions, Translation, add_point!,
                            check_point_ephemeris

using Ephemerides: EphemerisProvider, 
                   ephem_vector3, ephem_vector6, ephem_vector9, ephem_vector12


"""
    add_point_ephemeris!(fr::FrameSystem{O, N}, eph::EphemerisProvider, 
        name::Symbol, id::Int) where {O, N}
    
Add a point from `Ephemerides.jl` provider.
"""
function add_point_ephemeris!(
    fr::FrameSystem{O, N}, eph::EphemerisProvider, name::Symbol, id::Int
) where {O, N}

    pid, axid = check_point_ephemeris(fr, eph, id)

    funs = FramePointFunctions{O, N}(
        t -> Translation{O}( ephem_vector3(eph, pid, id, t) ),
        t -> Translation{O}( ephem_vector6(eph, pid, id, t) ),
        t -> Translation{O}( ephem_vector9(eph, pid, id, t) ),
        t -> Translation{O}( ephem_vector12(eph, pid, id, t) )
    )
    return add_point!(fr, name, id, axid, funs, pid)    
end

end