module EphemeridesExt 

import FrameTransformations: add_point_ephemeris!
using FrameTransformations: FrameSystem, SVectorNT, add_point!, FramePointFunctions,
                            check_retrieve, POINT_CLASSID_DYNAMIC

using Ephemerides: EphemerisProvider, ephem_vector3, ephem_vector6, ephem_vector9, ephem_vector12

function add_point_ephemeris!(
    frames::FrameSystem{O, N}, eph::EphemerisProvider, name::Symbol, id::Int
) where {O, N}

    parentid, axesid = check_retrieve(frames, eph, id)

    funs = FramePointFunctions{O, N}(
        t -> SVectorNT{3O, N}( ephem_vector3(eph, parentid, id, t) ),
        t -> SVectorNT{3O, N}( ephem_vector6(eph, parentid, id, t) ),
        t -> SVectorNT{3O, N}( ephem_vector9(eph, parentid, id, t) ),
        t -> SVectorNT{3O, N}( ephem_vector12(eph, parentid, id, t) )
    )

    add_point!(frames, name, id, axesid, POINT_CLASSID_DYNAMIC, funs, parentid)    
    
end

end