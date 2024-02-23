
export AXESID_ICRF, AXESID_GCRF, 
       AXESID_ECL2000, AXESID_EME2000, 
       AXESID_MOONME_DE421, AXESID_MOONPA_DE421, AXESID_MOONPA_DE440


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
    AXESID_MEME2000 
   
NAIF Axes ID for the Mean Ecliptic Equinox of J2000 (ECL2000).
"""
const AXESID_ECL2000 = 17

"""
    AXESID_EME2000

Axes ID for the Mean Dynamical Equator and Equinox of J2000.0. 

!!! note 
    In SPICE the J2000 and ICRF axes are considered equal, thus there exist no 
    specific NAIF ID for the EME2000 axes. 22 has been chosen because it is the 
    first unassigned axes ID among the built-in SPICE frames. 
"""
const AXESID_EME2000 = 22

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