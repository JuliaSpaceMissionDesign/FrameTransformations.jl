
@interface function add_point_ephemeris!(
    ::FrameSystem{O, N}, ::AbstractEphemerisProvider, ::Symbol, ::Int
) where {O, N} end

@interface function add_axes_ephemeris!(
    ::FrameSystem{O, N}, ::AbstractEphemerisProvider, ::Symbol, ::Int, ::Symbol
) where {O, N} end

function add_point_ephemeris!(
    frames::FrameSystem{O, N}, eph::AbstractEphemerisProvider, book::Dict{Int, Symbol}
) where {O, N}
    records = ephem_position_records(eph) 
    for id in sort(ephem_available_points(eph))
        !haskey(book, id) && throw(KeyError("Cannot find point with ID $id in the names book"))

        if id == 0 
            if length(get_points(frames).nodes) == 0 
                # No point registered in the frame system 
                add_point_root!(frames, book[id], id, AXESID_ICRF)
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

function check_point_ephemeris(
    frames::FrameSystem{O, N}, eph::AbstractEphemerisProvider, id::Int
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

    return parentid, axesid
end

function check_axes_ephemeris(
    frames::FrameSystem{O, N}, eph::AbstractEphemerisProvider, id::Int
) where {O, N}

    # Check that the kernels contain orientation data for the given axes ID 
    if !(id in ephem_available_axes(eph))
        throw(
            ArgumentError(
                "Orientation data for AXESID $id is not available " *
                "in the kernels loaded in the given FrameSystem.",
            ),
        )
    end

    # Retrieve the parent from the ephemeris data 
    orient_records = ephem_orient_records(eph)
    parentid = nothing

    for or in orient_records
        if or.target == id
            if isnothing(parentid)
                parentid = or.axes
            elseif parentid != or.axes
                throw(
                    ErrorException(
                        "UnambiguityError: at least two set of orientation data " *
                        "with different centers are available for axes with ID $id.",
                    ),
                )
            end
        end
    end

    # Check that the default parent is available in the FrameSystem! 
    if !has_axes(frames, parentid)
        throw(
            ArgumentError(
                "Orientation data for axes with ID $id is available " *
                "with respect to axes with ID $parentid, which has not yet been defined " *
                "in the given FrameSystem.",
            ),
        )
    end

    return parentid
end

# Functions for an easier definition\handling of ephemeris axes
@inline function triplet_to_rot3(θ, seq::Symbol)
    @inbounds angle_to_dcm(θ[1], θ[2], θ[3], seq)
end

@inline function triplet_to_rot6(θ, seq::Symbol)
    return _3angles_to_rot3(θ, seq), _3angles_to_δdcm(θ, seq)
end

@inline function triplet_to_rot9(θ, seq::Symbol)
    return (_3angles_to_rot3(θ, seq), _3angles_to_δdcm(θ, seq), _3angles_to_δ²dcm(θ, seq))
end

@inline function triplet_to_rot12(θ, seq::Symbol)
    return (
        _3angles_to_rot3(θ, seq), _3angles_to_δdcm(θ, seq), 
        _3angles_to_δ²dcm(θ, seq), _3angles_to_δ³dcm(θ, seq),
    )
end