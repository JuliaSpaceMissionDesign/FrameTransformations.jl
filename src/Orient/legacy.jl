
"""
    AXESID_B1950

NAIF Axes ID for the Mean Equator and Dynamical Equinox of the Besselian year 1950. 
"""
const AXESID_B1950 = 2

"""
    AXESID_FK4

NAIF Axes ID for the FK4 reference frame.
"""
const AXESID_FK4 = 3

"""
    AXESID_GALACTIC

NAIF Axes ID for the Galactic System II reference frame.
"""
const AXESID_GALACTIC = 13

"""
    AXESID_ECLIPB1950

NAIF Axes ID for the Mean Ecliptic Equinox of B1950.
"""
const AXESID_ECLIPB1950 = 18

"""
    DCM_MEME2000_TO_B1950

DCM for the rotation from the Mean Equator and Equinox of J2000.0 (`MEME2000`)
to the Mean Equator and Dynamical Equinox of B1950 (`B1950`).

!!! note 
    This rotation is obtained by precessing the J2000 frame backwards from Julian year 2000 
    to Besselian year 1950, using the 1976 IAU precession model. The rotation values 
    are taken from the SPICE toolkit. 

### References 
- SPICE [Library](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/chgirf.html)
"""
const DCM_MEME2000_TO_B1950 = angle_to_dcm(
    deg2rad(1153.04066200330/3600), -deg2rad(1002.26108439117/3600), 
    deg2rad(1152.84248596724/3600), :ZYZ
)

"""
    DCM_B1950_TO_ECLIPB1950

DCM for the rotation from the Mean Equator and Dynamical Equinox of B1950 (`MEMEB1950`) to
the Mean Ecliptic Equinox of B1950. This corresponds to the transformation `B1950 -> ECLIPB1950` 
in the SPICE toolkit, and uses the mean obliquity of the ecliptic from the IAU 1976 theory.

### References 
- SPICE [Library](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/chgirf.html)
"""
const DCM_B1950_TO_ECLIPB1950 = angle_to_dcm(deg2rad(84404.836/3600), :X)

"""
    DCM_B1950_TO_FK4

DCM for the rotation from the Mean Equator and Dynamical Equinox of B1950 (`MEMEB1950`) ot 
the FK4 reference frame. 

!!! note 
    The FK4 reference frame is obtained from the B1950 frame by applying the equinox offset 
    determined by Fricke.

### References 
- SPICE [Library](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/chgirf.html)
"""
const DCM_B1950_TO_FK4 = angle_to_dcm(deg2rad(0.525/3600), :Z)

"""
    DCM_FK4_TO_GALACTIC

DCM for the rotation from the FK4 frame to the Galactic System II reference frame.

!!! note 
    As the SPICE toolkit, we assume that this rotation is derived from the FK4 frame.

### References 
- SPICE [Library](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/spicelib/chgirf.html)
"""
const DCM_FK4_TO_GALACTIC = angle_to_dcm(deg2rad(282.25), deg2rad(62.6), deg2rad(327), :ZXZ)