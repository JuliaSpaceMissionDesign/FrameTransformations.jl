
export DCM_ICRF_TO_J2000_BIAS, 
       DCM_ICRF_TO_ECLIPJ2000,
       DCM_J2000_TO_ECLIPJ2000


""" 
    AXESID_ICRF 
   
NAIF Axes ID for the International Celestial Reference Frame (ICRF)
"""
AXESID_ICRF = 1

"""
    AXESID_MEME2000

Axes ID for the Mean Dynamical Equator and Equinox of J2000.0. 

!!! note 
    In SPICE the J2000 and ICRF axes are considered equal, thus there exist no 
    specific NAIF ID for the MEME2000 axes. 22 has been chosen because it is the 
    first unassigned axes ID among the built-in SPICE frames. 
"""
AXESID_MEME2000 = 22

""" 
    AXESID_ECLIPJ2000 
   
NAIF Axes ID for the Mean Ecliptic Equinox of J2000 (ECLIPJ2000) 
"""
AXESID_ECLIPJ2000 = 17


# --------------------------------------------------------
# DCMs
# --------------------------------------------------------


"""
    DCM_ICRF_TO_J2000_BIAS

DCM for the rotation from the International Celestial Reference Frame (`ICRF`) and the 
Mean Dynamical Equator and Equinox of J2000.0 (`MEME2000`).

### References
- Hilton, James L., and Catherine Y. Hohenkerk. -- Rotation matrix from the mean 
dynamical equator and equinox at J2000. 0 to the ICRS. -- Astronomy & Astrophysics 
413.2 (2004): 765-770. DOI: [10.1051/0004-6361:20031552](https://www.aanda.org/articles/aa/pdf/2004/02/aa3851.pdf)
- [SOFA docs](https://www.iausofa.org/2021_0512_C/sofa/sofa_pn_c.pdf)
"""
const DCM_ICRF_TO_J2000_BIAS = orient_precession_bias(iau2006a, 0.0)


"""
    DCM_J2000_TO_ECLIPJ2000

DCM for the rotation from the Mean Dynamical Equator of J2000 (`MEME2000`) to the 
Mean Ecliptic Equinox. This corresponds to the transformation `J2000 -> ECLIPJ2000` 
in the SPICE toolkit.
"""
const DCM_J2000_TO_ECLIPJ2000 = angle_to_dcm(orient_obliquity(iau2006a, 0.0), :X)


"""
    DCM_ICRF_TO_ECLIPJ2000

DCM for the rotation from the International Celestial Reference Frame (`ICRF`) to the 
Mean Ecliptic Equinox of J2000 (`ECLIPJ2000`).
"""
const DCM_ICRF_TO_ECLIPJ2000 = DCM_ICRF_TO_J2000_BIAS * DCM_J2000_TO_ECLIPJ2000


# --------------------------------------------------------
# TRANSFORMATIONS
# --------------------------------------------------------


"""
    orient_icrf_to_mememod(t::Number)

Compute the rotation matrix from the International Celestial Reference Frame (ICRF) to 
the Mean Equinox Mean Equator of Date at time `t`, expressed in TT seconds since [`J2000`](@ref).
"""
function orient_icrf_to_mememod(t::Number)
    # convert TT seconds since J2000 to TT centuries since J2000
    T = t/Tempo.CENTURY2SEC

    # fw_angles holds independent on the IAU Model! 
    γ, ϕ, ψ, ε = fw_angles(iau2006b, T)
    R = fw_matrix(γ, ϕ, ψ, ε)
    return R
end
