export add_axes_pa440!, add_axes_pa421!, add_axes_me421!

"""
	add_axes_pa440!(frames, axes::AbstractFrameAxes) 

Add `axes` as a set of ephemeris axes representing the DE440 Moon's Principal Axes (PA) to 
`frames`. The libration angles are extracted from the ephemeris kernels loaded within `frames`, 
an error is thrown if such orientation data is not available. 

!!! warning 
    The parent axes are automatically set to the ICRF (ID = $(Orient.AXESID_ICRF)). If such 
    axes are not registered in the frame system, an error is thrown.

!!! warning
	To properly read the ephemeris kernels, the ID associated to the input `axes` must 
	NAIF's FRAME ID for the Moon PA DE440 axes ($(Orient.AXESID_MOONPA_DE440)).
    
----

    add_axes_pa440!(frames, name::Symbol, axesid::Int = Orient.AXESID_MOONPA_DE440)

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) type 
via the [`@axes`](@ref) macro.

### See also 
See also [`Orient.AXESID_MOONPA_DE440`](@ref), [`Orient.orient_rot3_icrf_to_pa440`](@ref)
and [`add_axes_pa421!`](@ref).
"""
@inline function add_axes_pa440!(frames::FrameSystem, axes::AbstractFrameAxes)
    return add_axes_pa440!(frames, axes_name(axes), axes_id(axes))
end

# Low-level function
function add_axes_pa440!(
    frames::FrameSystem, name::Symbol, axesid::Int=Orient.AXESID_MOONPA_DE440
)

    if !(has_axes(frames, Orient.AXESID_ICRF))
        throw(
            ErrorException(
                "The DE440 Moon Principal Axes (PA) can only be defined w.r.t. the" * 
                " ICRF (ID = $(Orient.AXESID_ICRF)), which is not defined" *
                " in the current frames graph."
            )
        )
    end

    # Throw a warning if the ID does not match the standard one.
    if axesid != Orient.AXESID_MOONPA_DE440
        throw(
            ArgumentError(
                "$name is aliasing an ID that is not the standard PA440" * 
                " ID ($(Orient.AXESID_MOONPA_DE440))."
            )
        )
    end

    return add_axes_ephemeris!(frames, name, axesid, :ZXZ)
end


"""
	add_axes_pa421!(frames, axes::AbstractFrameAxes) 

Add `axes` as a set of ephemeris axes representing the DE421 Moon's Principal Axes (PA) to 
`frames`. The libration angles are extracted from the ephemeris kernels loaded within `frames`, 
an error is thrown if such orientation data is not available. 

!!! warning 
    The parent axes are automatically set to the ICRF (ID = $(Orient.AXESID_ICRF)). If such 
    axes are not registered in the frame system, an error is thrown.

!!! warning
	To properly read the ephemeris kernels, the ID associated to the input `axes` must match 
	NAIF's FRAME ID for the Moon PA DE421 axes ($(Orient.AXESID_MOONPA_DE421)).

----

    add_axes_pa421!(frames, name::Symbol, axesid::Int = Orient.AXESID_MOONPA_DE421)

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) type 
via the [`@axes`](@ref) macro.

### See also 
See also [`Orient.AXESID_MOONPA_DE421`](@ref), [`Orient.orient_rot3_icrf_to_pa421`](@ref), 
[`add_axes_pa440!`](@ref), and [`add_axes_me421!`](@ref)
"""
@inline function add_axes_pa421!(frames::FrameSystem, axes::AbstractFrameAxes)
    return add_axes_pa421!(frames, axes_name(axes), axes_id(axes))
end

# Low-level function
function add_axes_pa421!(
    frames::FrameSystem, name::Symbol, axesid::Int=Orient.AXESID_MOONPA_DE421
)

    if !(has_axes(frames, Orient.AXESID_ICRF))
        throw(
            ErrorException(
                "The DE421 Moon Principal Axes (PA) can only be defined w.r.t. the" * 
                " ICRF (ID = $(Orient.AXESID_ICRF)), which is not defined" *
                " in the current frames graph."
            )
        )
    end

    # Throw a warning if the ID does not match the standard one.
    if axesid != Orient.AXESID_MOONPA_DE421
        throw(
            ArgumentError(
                "$name is aliasing an ID that is not the standard PA421" * 
                " ID ($(Orient.AXESID_MOONPA_DE421))."
            )
        )
    end

    return add_axes_ephemeris!(frames, name, axesid, :ZXZ)

end


"""
	add_axes_me421!(frames, axes::AbstractFrameAxes, parent) 

Add `axes` as fixed offset axes representing the DE421 Moon's Mean Earth/Mean Rotation (ME) 
to `frames`.

!!! warning 
    The `parent` set of axes must either the DE440 Principal Axes (PA440, ID = 
    $(Orient.AXESID_MOONPA_DE440)) or the DE421 Principal Axes (PA421, ID = 
    $(Orient.AXESID_MOONPA_DE421)), otherwise an error is thrown. Depending on that, the 
    relative axes orientation will be automatically selected by this function. 

----

    add_axes_me421!(frames, name::Symbol, parentid::Int, axesid::Int = Orient.AXESID_MOONME_DE421)

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) type 
via the [`@axes`](@ref) macro.

### See also 
See also [`add_axes_pa440!`](@ref), and [`add_axes_pa421!`](@ref), 
[`Orient.DCM_MOON_PA421_TO_ME421`](@ref) and [`Orient.DCM_MOON_PA421_TO_ME421`](@ref), 
"""
@inline function add_axes_me421!(frames::FrameSystem, axes::AbstractFrameAxes, parent)
    return add_axes_me421!(frames, axes_name(axes), axes_alias(parent), axes_id(axes))
end

# Low-level function
function add_axes_me421!(
    frames::FrameSystem, name::Symbol, parentid::Int, axesid::Int=Orient.AXESID_MOONME_DE421
)

    if parentid == Orient.AXESID_MOONPA_DE421
        dcm = Orient.DCM_MOON_PA421_TO_ME421

    elseif parentid == Orient.AXESID_MOONPA_DE440
        dcm = Orient.DCM_MOON_PA440_TO_ME421

    else

        throw(
            ArgumentError(
                "The DE421 Mean Earth/Mean Rotation (ME) axes cannot be defined w.r.t." *
                " axes with ID $parentid. Only the DE440 (ID = $(Orient.AXESID_MOONPA_DE440))" * 
                " or DE421 (ID = $(Orient.AXESID_MOONPA_DE421)) Moon Principal Axes are" *
                " accepted as parent axes.",
            ),
        )
    end

    if axesid != Orient.AXESID_MOONME_DE421
        @warn "$(name) is aliasing an ID that is not the standard ME421" * 
            "ID ($(Orient.AXESID_MOONME_DE421))."
    end

    return add_axes_fixedoffset!(frames, name, axesid, parentid, dcm)

end
