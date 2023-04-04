export DCM_MOON_PA440_TO_ME421,
    DCM_MOON_PA430_TO_ME430,
    DCM_MOON_PA430_TO_ME421,
    DCM_MOON_PA421_TO_ME421,
    orient_rot3_icrf_to_pa440,
    orient_rot3_icrf_to_pa421

"""
    AXESID_MOONPA_DE421 

NAIF axes id for the DE421 Moon Principal Axes (PA421).
"""
const AXESID_MOONPA_DE421 = 31006

""" 
    AXESID_MOONME_DE421
    
NAIF axes id for the DE421 Moon Mean Earth/Mean Rotation axes  (ME421).
"""
const AXESID_MOONME_DE421 = 31007

""" 
    AXESID_MOONPA_DE440 

NAIF Axes id for the DE440 Moon Principal Axes (PA440).
"""
const AXESID_MOONPA_DE440 = 31008

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
input time `t`, expressed in seconds since [`J2000`](@ref). 

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
    ephem_orient_order!(y, eph, DJ2000, t / Tempo.DAY2SEC, AXESID_MOONPA_DE440, 0)
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
input time `t`, expressed in seconds since [`J2000`](@ref). 

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
    ephem_orient_order!(y, eph, DJ2000, t / Tempo.DAY2SEC, AXESID_MOONPA_DE421, 0)
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
