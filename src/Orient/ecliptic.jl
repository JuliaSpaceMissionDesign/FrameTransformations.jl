
export DCM_ICRF_TO_ECLIPJ2000, DCM_MEME2000_TO_ECLIPJ2000

""" 
    AXESID_ICRF 
   
NAIF Axes ID for the International Celestial Reference Frame (ICRF)
"""
const AXESID_ICRF = 1

"""
    AXESID_GCRF 

Axes ID for the Geocentric Celestial Reference Frame (GCRFF)

!!! note 
    Although the ICRF and GCRF axes are identical, they are based upon a different 
    timescale. A different ID is here assigned to provide a robust way of distinguishing 
    between the two. 23 has been chosen because it is one the unassigned axes ID among the 
    built-in SPICE frames.
"""
const AXESID_GCRF = 23

""" 
    AXESID_ECLIPJ2000 
   
NAIF Axes ID for the Mean Ecliptic Equinox of J2000 (ECLIPJ2000) 
"""
const AXESID_ECLIPJ2000 = 17

# --------------------------------------------------------
# DCMs
# --------------------------------------------------------

"""
    DCM_MEME2000_TO_ECLIPJ2000

DCM for the rotation from the Mean Equator and Equinox of J2000 (`MEME2000`) to the 
Mean Ecliptic Equinox. This corresponds to the transformation `J2000 -> ECLIPJ2000` 
in the SPICE toolkit, and uses the mean obliquity of the ecliptic from the IAU 1976 theory.
"""
const DCM_MEME2000_TO_ECLIPJ2000 = angle_to_dcm(orient_obliquity(iau1980, 0.0), :X)

"""
    DCM_ICRF_TO_ECLIPJ2000

DCM for the rotation from the International Celestial Reference Frame (`ICRF`) to the 
Mean Ecliptic Equinox of J2000 (`ECLIPJ2000`).
"""
const DCM_ICRF_TO_ECLIPJ2000 = DCM_ICRF_TO_MEME2000 * DCM_MEME2000_TO_ECLIPJ2000
