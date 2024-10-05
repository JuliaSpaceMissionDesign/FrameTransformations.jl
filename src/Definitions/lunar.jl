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
    add_axes_pa440!(frames, eph::AbstractEphemerisProvider, name::Symbol, 
        id::Int=AXESID_MOONPA_DE440)

Add DE440 Moon's Principal Axes (PA) axes to `frames`. The libration angles are extracted 
from the `eph` ephemeris kernels, an error is thrown if such orientation data is not available. 

!!! warning 
    The parent axes are automatically set to the ICRF (ID = $(AXESID_ICRF)). If such 
    axes are not registered in the frame system, an error is thrown.

!!! warning
	To properly read the ephemeris kernels, the ID associated to the input `axes` must match 
	NAIF's FRAME ID for the Moon PA DE440 axes ($(AXESID_MOONPA_DE440)).
"""
function add_axes_pa440!(
    frames::FrameSystem, eph::AbstractEphemerisProvider, name::Symbol,
    id::Int=AXESID_MOONPA_DE440
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
    if id != AXESID_MOONPA_DE440
        throw(
            ArgumentError(
                "$name is aliasing an ID that is not the standard PA440" *
                " ID ($(AXESID_MOONPA_DE440))."
            )
        )
    end

    return add_axes_ephemeris!(frames, eph, name, id, :ZXZ)
end

"""
    add_axes_pa421!(frames, eph::AbstractEphemerisProvider, name::Symbol, 
        id::Int=AXESID_MOONPA_DE421)

Add DE421 Moon's Principal Axes (PA) axes to `frames`. The libration angles are extracted 
from the `eph` ephemeris kernels, an error is thrown if such orientation data is not available. 

!!! warning 
    The parent axes are automatically set to the ICRF (ID = $(AXESID_ICRF)). If such 
    axes are not registered in the frame system, an error is thrown.

!!! warning
	To properly read the ephemeris kernels, the ID associated to the input `axes` must match 
	NAIF's FRAME ID for the Moon PA DE421 axes ($(AXESID_MOONPA_DE421)).
"""
function add_axes_pa421!(
    frames::FrameSystem, eph::AbstractEphemerisProvider, name::Symbol,
    id::Int=AXESID_MOONPA_DE421
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
    if id != AXESID_MOONPA_DE421
        throw(
            ArgumentError(
                "$name is aliasing an ID that is not the standard PA421" *
                " ID ($(AXESID_MOONPA_DE421))."
            )
        )
    end

    return add_axes_ephemeris!(frames, eph, name, id, :ZXZ)
end

"""
    add_axes_me421!(frames, name::Symbol, parentid::Int, axesid::Int=AXESID_MOONME_DE421)

Add DE421 Moon's Mean Earth/Mean Rotation (ME) axes to `frames`.

!!! warning 
    The `parent` set of axes must either the DE440 Principal Axes (PA440, ID = 
    $(AXESID_MOONPA_DE440)) or the DE421 Principal Axes (PA421, ID = 
    $(AXESID_MOONPA_DE421)), otherwise an error is thrown. Depending on that, the 
    relative axes orientation will be automatically selected by this function. 
"""
function add_axes_me421!(
    frames::FrameSystem, name::Symbol, parent, id::Int=AXESID_MOONME_DE421
)
    pid = axes_id(frames, parent)
    if pid == AXESID_MOONPA_DE421
        dcm = DCM_MOON_PA421_TO_ME421

    elseif pid == AXESID_MOONPA_DE440
        dcm = DCM_MOON_PA440_TO_ME421

    else

        throw(
            ArgumentError(
                "The DE421 Mean Earth/Mean Rotation (ME) axes cannot be defined w.r.t." *
                " axes with ID $pid. Only the DE440 (ID = $(AXESID_MOONPA_DE440))" *
                " or DE421 (ID = $(AXESID_MOONPA_DE421)) Moon Principal Axes are" *
                " accepted as parent axes.",
            ),
        )
    end

    if id != AXESID_MOONME_DE421
        @warn "$(name) is aliasing an ID that is not the standard ME421" *
              "ID ($(AXESID_MOONME_DE421))."
    end

    return add_axes_fixedoffset!(frames, name, id, pid, dcm)
end
