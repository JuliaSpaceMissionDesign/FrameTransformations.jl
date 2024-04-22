module EphemeridesExt 

import FrameTransformations: add_point_ephemeris!
using FrameTransformations: FrameSystem, SVectorNT,
                            add_point!, get_points, FramePointFunctions,
                            add_point_root!

using JSMDInterfaces.Ephemeris: ephem_available_points, ephem_position_records
using Ephemerides: EphemerisProvider, 
                   ephem_vector3, ephem_vector6, ephem_vector9, ephem_vector12

function add_point_ephemeris!(
    frames::FrameSystem{O, N}, eph::EphemerisProvider, name::Symbol, id::Int
) where {O, N}

    # Check that the kernels contain the ephemeris data for the given naifid
    if !(id in ephem_available_points(eph))
        throw(
            ErrorException(
                "Ephemeris data for ID $naifid is not available in the kernels.",
            ),
        )
    end

    # Retrieve the ephemerides position records (i.e., the segment descriptors)
    pos_records = ephem_position_records(frames.eph)

    # Retrieve the parent point and the point axes from the ephemeris data 
    parentid, parent_found = 0, false
    axesid, axes_found = 0, false

    for pr in pos_records
        if pr.target == id
            
            # Update the parent point ID
            if !parent_found
                parentid = pr.center
                parent_found = true 

            elseif parentid != pr.center
                throw(
                    ErrorException(
                        "UnambiguityError: at least two set of data with different" * 
                        " centers are available for point with ID $id.",
                    ),
                )

            end
            
            # Update the reference axes ID
            if !axes_found
                axesid = pr.axes 
                axes_found = true 

            elseif axesid != pr.axes 
                throw(
                    ErrorException(
                        "UnambiguityError: at least two set of data with different axes" *
                        " are available for point with ID $id.",
                    ),
                )
            end

        end
    end

     # Check that the default parent is available in the FrameSystem
     if !has_point(frames, parentid)
        throw(
            ErrorException(
                "Ephemeris data for point with ID $id is available with respect" * 
                " to point with ID $parentid, which has not yet been defined" *
                " in the input frame system.",
            ),
        )
    end

    # The parent point is automatically inferred from the ephemeris kernels so it will 
    # always have available ephemeris data

    # Checks if the axes are known to the frame system. 
    # This check is also performed by add_point!, but it is reported here because 
    # it provides more specific information for ephemeris points 
    if !has_axes(frames, axesid)
        throw(
            ErrorException(
                "Ephemeris data for point with ID $naifid is expressed in a set" *
                " of axes with ID $axesid, which are yet to be defined in the" *
                " input frame system.",
            ),
        )
    end

    funs = FramePointFunctions{O, N}(
        t -> SVectorNT{3O, N}( ephem_vector3(eph, parentid, id, t) ),
        t -> SVectorNT{3O, N}( ephem_vector6(eph, parentid, id, t) ),
        t -> SVectorNT{3O, N}( ephem_vector9(eph, parentid, id, t) ),
        t -> ephem_vector12(eph, parentid, id, t) 
    )

    add_point!(frames, name, id, axesid, 3, funs, parentid) # TODO: change class to a const 
    nothing 

end

function add_point_ephemeris!(
    frames::FrameSystem{O, N}, eph::EphemerisProvider, book::Dict{Int, Symbol}
) where {O, N}
    records = ephem_position_records(eph) 
    for id in sort(ephem_available_points(eph))
        !haskey(book, id) && throw(KeyError("Cannot find point with ID $id in the names book"))

        if id == 0 
            if length(get_points(frames).nodes) == 0 
                # No point registered in the frame system 
                add_point_root!(book[id], id, 1) # TODO: change this 1 to AXESID_ICRF
            end
        else
            # Assume that the SSB is registered in the frame system
            it = findfirst(rec -> rec.target == id, records)
            rec = records[it]
            add_point_ephemeris!(frames, eph, book[id], id)
        end
    end
    nothing 
end

end