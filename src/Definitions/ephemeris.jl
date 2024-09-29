
"""
    add_point_ephemeris!(fr::FrameSystem{O, T}, eph::AbstractEphemerisProvider, 
        name::Symbol, id::Int) where {O, T} 

Add a point coming from an `AbstractEphemerisProvider` subtype. 
The point is identifies by the `id` and a have a user defined `name`.

!!! warning 
    This is an interface only, concrete subtypes of `AbstractEphemerisProvider` requires 
    an proper implementation.
"""
@interface function add_point_ephemeris!(
    ::FrameSystem{O,T}, ::AbstractEphemerisProvider, ::Symbol, ::Int
) where {O,T} end

"""
    add_point_ephemeris!(fr, eph::AbstractEphemerisProvider, book::Dict{Int, Symbol})

Add all points found in the `eph` and using id-names relationships specified in `book`.
"""
function add_point_ephemeris!(
    fr::FrameSystem{O,T}, eph::AbstractEphemerisProvider, book::Dict{Int,Symbol}
) where {O,T}
    records = ephem_position_records(eph)
    for id in sort(ephem_available_points(eph))
        !haskey(book, id) && @warn "Cannot find point with ID $id in the names book" continue

        if id == 0
            if length(points_graph(fr).nodes) == 0
                # No point registered in the frame system 
                add_point!(fr, book[id], id, AXESID_ICRF)
            end
        else
            # Assume that the SSB is registered in the frame system
            it = findfirst(rec -> rec.target == id, records)
            rec = records[it]
            add_point_ephemeris!(fr, eph, book[id], id)
        end
    end
    nothing
end

function check_point_ephemeris(
    fr::FrameSystem{O,T}, eph::AbstractEphemerisProvider, id::Int
) where {O,T}

    # Check that the kernels contain the ephemeris data for the given naifid
    if !(id in ephem_available_points(eph))
        throw(
            ErrorException(
                "Ephemeris data for ID $naifid is not available in the kernels.",
            ),
        )
    end

    # Retrieve the ephemerides position records (i.e., the segment descriptors)
    pos_records = ephem_position_records(eph)

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
    if !has_point(fr, parentid)
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
    if !has_axes(fr, axesid)
        throw(
            ErrorException(
                "Ephemeris data for point with ID $naifid is expressed in a set" *
                " of axes with ID $axesid, which are yet to be defined in the" *
                " input frame system.",
            ),
        )
    end

    return parentid, axesid
end
