export DCM_MOONPA440_TO_MER421, 
       DCM_MOONPA430_TO_MER430, 
       DCM_MOONPA430_TO_MER421,
       DCM_MOONPA421_TO_MER421,
       orient_rot3_icrf_to_pa440, 
       orient_rot3_icrf_to_pa421

""" 
    AXESID_MOONPA_DE440 

NAIF Axes id for the DE440 Moon Principal Axes.
"""
AXESID_MOONPA_DE440 = 31008


"""
    AXESID_MOONPA_DE421 

NAIF axes id for the DE421 Moon Principal Axes 
"""
AXESID_MOONPA_DE421 = 31006

"""
    DCM_MOONPA440_TO_MER421

DCM for the rotation from the Moon Principal Axis 440 (PA440) to the Moon Mean 
Earth/Mean Rotation (MER) axes.

### References 
 - Park, S. R. et al. (2021), _The JPL Planetary and Lunar Ephemerides DE440 and DE441_,
 [DOI: 10.3847/1538-3881/abd414](https://doi.org/10.3847/1538-3881/abd414) 

"""
DCM_MOONPA440_TO_MER421 = angle_to_dcm(arcsec2rad(-67.8526), arcsec2rad(-78.6944), 
                                       arcsec2rad(-0.2785), :ZYX)


""" 
    DCM_MOONPA430_TO_MER430

DCM for the rotation from the Moon Principal Axis 430 (PA430) to the Moon Mean 
Earth/Mean Rotation (MER430) axes.

### References 
- Folkner M. William et al. (2014), _The Planetary and Lunar EphemeridesDE430 and DE431_

- J. G. Williams et al. (2013), _DE430 Lunar Orbit, Physical Librations, and Surface Coordinates_,
[DE430 Lunar Ephemeris and Orientation](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430_moon_coord.pdf) 

"""
DCM_MOONPA430_TO_MER430 = angle_to_dcm(arcsec2rad(-67.573), arcsec2rad(-78.58), 
                                       arcsec2rad(-0.285), :ZYX)


""" 
    DCM_MOONPA430_TO_MER421

DCM for the rotation from the Moon Principal Axis 430 (PA430) to the Moon Mean 
Earth/Mean Rotation (MER421) axes.

### References 
- Folkner M. William et al. (2014), _The Planetary and Lunar EphemeridesDE430 and DE431_

- J. G. Williams et al. (2013), _DE430 Lunar Orbit, Physical Librations, and Surface Coordinates_,
[DE430 Lunar Ephemeris and Orientation](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430_moon_coord.pdf) 

"""
DCM_MOONPA430_TO_MER421 = angle_to_dcm(arcsec2rad(-67.737), arcsec2rad(-78.627), 
                                       arcsec2rad(-0.295), :ZYX)


"""
    DCM_MOONPA421_TO_MER421

DCM for the rotation from the Moon Principal Axis 421 (PA421) to the Moon Mean 
Earth/Mean Rotation (MER421) axes.

### References 
- J. G. Williams et al. (2008), _DE421 Lunar Orbit, Physical Librations, and Surface Coordinates_,
[DE421 Lunar Ephemeris and Orientation](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/de421_lunar_ephemeris_and_orientation.pdf) 

"""
DCM_MOONPA421_TO_MER421 = angle_to_dcm(arcsec2rad(-67.92), arcsec2rad(-78.56),
                                       arcsec2rad(-0.3), :ZYX)

                                       
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
        throw(ErrorException(
            "Orientation data for the DE440 Moon Principal Axes is not available in the kernels "*
            "provided."
        ))
    end 

    y = @MVector zeros(3)
    ephem_orient_order!(y, eph, DJ2000, t/Tempo.DAY2SEC, AXESID_MOONPA_DE440, 0)
    angle_to_dcm(y[1], y[2], y[3], :ZXZ)
end

"""
    orient_rot3_icrf_to_pa440(eph::AbstractEphemerisProvider, ep::Epoch)

Compute the rotation matrix from the ICRF to the DE440 Moon's Principal Axes at the input 
epoch `ep`.
"""
function orient_rot3_icrf_to_pa440(eph::AbstractEphemerisProvider, ep::Epoch)
    orient_rot3_icrf_to_pa440(eph, j2000s(convert(TDB, ep)))
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
        throw(ErrorException(
            "Orientation data for the DE421 Moon Principal Axes is not available "*
            "in the kernels provided."
        ))
    end 

    y = @MVector zeros(3)
    ephem_orient_order!(y, eph, DJ2000, t/Tempo.DAY2SEC, AXESID_MOONPA_DE421, 0)
    angle_to_dcm(y[1], y[2], y[3], :ZXZ)
end


"""
    orient_rot3_icrf_to_pa421(eph::AbstractEphemerisProvider, ep::Epoch)

Compute the rotation matrix from the ICRF to the DE421 Moon's Principal Axes at the input 
epoch `ep`.
"""
function orient_rot3_icrf_to_pa421(eph::AbstractEphemerisProvider, ep::Epoch)
    orient_rot3_icrf_to_pa421(eph, j2000s(convert(TDB, ep)))
end
