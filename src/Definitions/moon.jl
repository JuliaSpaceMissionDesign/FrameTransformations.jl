export  DCM_MOON_PA440_TO_ME421, 
        DCM_MOON_PA430_TO_ME430, 
        DCM_MOON_PA430_TO_ME421, 
        DCM_MOON_PA421_TO_ME421,
        add_axes_pa440!, 
        add_axes_pa421!, 
        add_axes_me421!

"""
    DCM_MOON_PA440_TO_ME421

DCM for the rotation from the Moon Principal Axis 440 (PA440) to the Moon Mean 
Earth/Mean Rotation DE421 (ME421) axes.

### References 
 - Park, S. R. et al. (2021), _The JPL Planetary and Lunar Ephemerides DE440 and DE441_,
    [DOI: 10.3847/1538-3881/abd414](https://doi.org/10.3847/1538-3881/abd414) 

"""
const DCM_MOON_PA440_TO_ME421 = angle_to_dcm(
    arcsec2rad(-67.8526), arcsec2rad(-78.6944), arcsec2rad(-0.2785), :ZYX
)

""" 
    DCM_MOON_PA430_TO_ME430

DCM for the rotation from the Moon Principal Axis 430 (PA430) to the Moon Mean 
Earth/Mean Rotation DE430 (ME430) axes.

### References 
- Folkner M. William et al. (2014), _The Planetary and Lunar EphemeridesDE430 and DE431_

- J. G. Williams et al. (2013), _DE430 Lunar Orbit, Physical Librations, and Surface Coordinates_,
    [DE430 Lunar Ephemeris and Orientation](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430_moon_coord.pdf) 

"""
const DCM_MOON_PA430_TO_ME430 = angle_to_dcm(
    arcsec2rad(-67.573), arcsec2rad(-78.58), arcsec2rad(-0.285), :ZYX
)

""" 
    DCM_MOON_PA430_TO_ME421

DCM for the rotation from the Moon Principal Axis 430 (PA430) to the Moon Mean 
Earth/Mean Rotation DE421 (ME421) axes.

### References 
- Folkner M. William et al. (2014), _The Planetary and Lunar EphemeridesDE430 and DE431_

- J. G. Williams et al. (2013), _DE430 Lunar Orbit, Physical Librations, and Surface Coordinates_,
    [DE430 Lunar Ephemeris and Orientation](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430_moon_coord.pdf) 

"""
const DCM_MOON_PA430_TO_ME421 = angle_to_dcm(
    arcsec2rad(-67.737), arcsec2rad(-78.627), arcsec2rad(-0.295), :ZYX
)

"""
    DCM_MOON_PA421_TO_ME421

DCM for the rotation from the Moon Principal Axis 421 (PA421) to the Moon Mean 
Earth/Mean Rotation DE421 (ME421) axes.

### References 
- J. G. Williams et al. (2008), _DE421 Lunar Orbit, Physical Librations, and Surface Coordinates_,
    [DE421 Lunar Ephemeris and Orientation](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/de421_lunar_ephemeris_and_orientation.pdf) 

"""
const DCM_MOON_PA421_TO_ME421 = angle_to_dcm(
    arcsec2rad(-67.92), arcsec2rad(-78.56), arcsec2rad(-0.3), :ZYX
)

""" 
    orient_rot3_icrf_to_pa440(eph::AbstractEphemerisProvider, t::Number)

Compute the rotation matrix from the ICRF to the DE440 Moon's Principal Axes at the given 
input time `t`, expressed in seconds since `J2000`. 

!!! warning 
    This function is not optimised for performance (it allocates!). The user is suggested 
    to retrieve the Principal axes orientation using the dedicated Frame System functions.

"""
function orient_rot3_icrf_to_pa440(eph::AbstractEphemerisProvider, t::Number)
    if !(AXESID_MOONPA_DE440 in ephem_available_axes(eph))
        throw(
            ErrorException(
                "Orientation data for the DE440 Moon Principal Axes is not available in the kernels " *
                "provided.",
            ),
        )
    end

    y = @MVector zeros(3)
    ephem_orient!(y, eph, DJ2000, t / Tempo.DAY2SEC, AXESID_MOONPA_DE440, AXESID_ICRF, 0)
    return angle_to_dcm(y[1], y[2], y[3], :ZXZ)
end

"""
    orient_rot3_icrf_to_pa440(eph::AbstractEphemerisProvider, ep::Epoch)

Compute the rotation matrix from the ICRF to the DE440 Moon's Principal Axes at the input 
epoch `ep`.
"""
function orient_rot3_icrf_to_pa440(eph::AbstractEphemerisProvider, ep::Epoch)
    return orient_rot3_icrf_to_pa440(eph, j2000s(convert(TDB, ep)))
end

""" 
    orient_rot3_icrf_to_pa421(eph::AbstractEphemerisProvider, t::Number)

Compute the rotation matrix from the ICRF to the DE421 Moon's Principal Axes at the given 
input time `t`, expressed in seconds since `J2000`. 

!!! warning 
    This function is not optimised for performance (it allocates!). The user is suggested 
    to retrieve the Principal axes orientation using the dedicated Frame System functions.

"""
function orient_rot3_icrf_to_pa421(eph::AbstractEphemerisProvider, t::Number)
    if !(AXESID_MOONPA_DE421 in ephem_available_axes(eph))
        throw(
            ErrorException(
                "Orientation data for the DE421 Moon Principal Axes is not available " *
                "in the kernels provided.",
            ),
        )
    end

    y = @MVector zeros(3)
    ephem_orient!(y, eph, DJ2000, t / Tempo.DAY2SEC, AXESID_MOONPA_DE421, AXESID_ICRF, 0)
    return angle_to_dcm(y[1], y[2], y[3], :ZXZ)
end

"""
    orient_rot3_icrf_to_pa421(eph::AbstractEphemerisProvider, ep::Epoch)

Compute the rotation matrix from the ICRF to the DE421 Moon's Principal Axes at the input 
epoch `ep`.
"""
function orient_rot3_icrf_to_pa421(eph::AbstractEphemerisProvider, ep::Epoch)
    return orient_rot3_icrf_to_pa421(eph, j2000s(convert(TDB, ep)))
end

"""
	add_axes_pa440!(frames, axes::AbstractFrameAxes) 

Add `axes` as a set of ephemeris axes representing the DE440 Moon's Principal Axes (PA) to 
`frames`. The libration angles are extracted from the ephemeris kernels loaded within `frames`, 
an error is thrown if such orientation data is not available. 

!!! warning 
    The parent axes are automatically set to the ICRF (ID = $(AXESID_ICRF)). If such 
    axes are not registered in the frame system, an error is thrown.

!!! warning
	To properly read the ephemeris kernels, the ID associated to the input `axes` must match 
	NAIF's FRAME ID for the Moon PA DE440 axes ($(AXESID_MOONPA_DE440)).
    
----

    add_axes_pa440!(frames, name::Symbol, axesid::Int = AXESID_MOONPA_DE440)

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) type 
via the [`@axes`](@ref) macro.

### See also 
See also [`AXESID_MOONPA_DE440`](@ref), [`orient_rot3_icrf_to_pa440`](@ref)
and [`add_axes_pa421!`](@ref).
"""
@inline function add_axes_pa440!(frames::FrameSystem, axes::AbstractFrameAxes)
    return add_axes_pa440!(frames, axes_name(axes), axes_id(axes))
end

# Low-level function
function add_axes_pa440!(
    frames::FrameSystem, name::Symbol, axesid::Int=AXESID_MOONPA_DE440
)
    if !(has_axes(frames, AXESID_ICRF))
        throw(
            ErrorException(
                "The DE440 Moon Principal Axes (PA) can only be defined w.r.t. the" * 
                " ICRF (ID = $(AXESID_ICRF)), which is not defined" *
                " in the current frames graph."
            )
        )
    end

    # Throw a warning if the ID does not match the standard one.
    if axesid != AXESID_MOONPA_DE440
        throw(
            ArgumentError(
                "$name is aliasing an ID that is not the standard PA440" * 
                " ID ($(AXESID_MOONPA_DE440))."
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
    The parent axes are automatically set to the ICRF (ID = $(AXESID_ICRF)). If such 
    axes are not registered in the frame system, an error is thrown.

!!! warning
	To properly read the ephemeris kernels, the ID associated to the input `axes` must match 
	NAIF's FRAME ID for the Moon PA DE421 axes ($(AXESID_MOONPA_DE421)).

----

    add_axes_pa421!(frames, name::Symbol, axesid::Int = AXESID_MOONPA_DE421)

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) type 
via the [`@axes`](@ref) macro.

### See also 
See also [`AXESID_MOONPA_DE421`](@ref), [`orient_rot3_icrf_to_pa421`](@ref), 
[`add_axes_pa440!`](@ref), and [`add_axes_me421!`](@ref)
"""
@inline function add_axes_pa421!(frames::FrameSystem, axes::AbstractFrameAxes)
    return add_axes_pa421!(frames, axes_name(axes), axes_id(axes))
end

# Low-level function
function add_axes_pa421!(
    frames::FrameSystem, name::Symbol, axesid::Int=AXESID_MOONPA_DE421
)
    if !(has_axes(frames, AXESID_ICRF))
        throw(
            ErrorException(
                "The DE421 Moon Principal Axes (PA) can only be defined w.r.t. the" * 
                " ICRF (ID = $(AXESID_ICRF)), which is not defined" *
                " in the current frames graph."
            )
        )
    end

    # Throw a warning if the ID does not match the standard one.
    if axesid != AXESID_MOONPA_DE421
        throw(
            ArgumentError(
                "$name is aliasing an ID that is not the standard PA421" * 
                " ID ($(AXESID_MOONPA_DE421))."
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
    $(AXESID_MOONPA_DE440)) or the DE421 Principal Axes (PA421, ID = 
    $(AXESID_MOONPA_DE421)), otherwise an error is thrown. Depending on that, the 
    relative axes orientation will be automatically selected by this function. 

----

    add_axes_me421!(frames, name::Symbol, parentid::Int, axesid::Int = AXESID_MOONME_DE421)

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) type 
via the [`@axes`](@ref) macro.

### See also 
See also [`add_axes_pa440!`](@ref), and [`add_axes_pa421!`](@ref), 
[`DCM_MOON_PA421_TO_ME421`](@ref) and [`DCM_MOON_PA421_TO_ME421`](@ref), 
"""
@inline function add_axes_me421!(frames::FrameSystem, axes::AbstractFrameAxes, parent)
    return add_axes_me421!(frames, axes_name(axes), axes_alias(parent), axes_id(axes))
end

# Low-level function
function add_axes_me421!(
    frames::FrameSystem, name::Symbol, parentid::Int, axesid::Int=AXESID_MOONME_DE421
)
    if parentid == AXESID_MOONPA_DE421
        dcm = DCM_MOON_PA421_TO_ME421

    elseif parentid == AXESID_MOONPA_DE440
        dcm = DCM_MOON_PA440_TO_ME421

    else

        throw(
            ArgumentError(
                "The DE421 Mean Earth/Mean Rotation (ME) axes cannot be defined w.r.t." *
                " axes with ID $parentid. Only the DE440 (ID = $(AXESID_MOONPA_DE440))" * 
                " or DE421 (ID = $(AXESID_MOONPA_DE421)) Moon Principal Axes are" *
                " accepted as parent axes.",
            ),
        )
    end

    if axesid != AXESID_MOONME_DE421
        @warn "$(name) is aliasing an ID that is not the standard ME421" * 
            "ID ($(AXESID_MOONME_DE421))."
    end

    return add_axes_fixedoffset!(frames, name, axesid, parentid, dcm)
end