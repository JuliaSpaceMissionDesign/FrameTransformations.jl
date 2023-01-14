
export DCM_ICRF2J2000_BIAS, 
       DCM_J20002ECLIPJ2000, 
       DCM_ICRF2ECLIPJ2000

# --------------------------------------------------------
# DCMs
# --------------------------------------------------------

"""
    DCM_ICRF2J2000_BIAS

DCM for the rotation from the International Celestial Reference Frame (`ICRF`) and the 
Mean Dynamical Equator and Equinox of J2000.0 (`MEME2000`).

### References
- Hilton, James L., and Catherine Y. Hohenkerk. -- Rotation matrix from the mean 
    dynamical equator and equinox at J2000. 0 to the ICRS. -- Astronomy & Astrophysics 
    413.2 (2004): 765-770. DOI: [10.1051/0004-6361:20031552](https://www.aanda.org/articles/aa/pdf/2004/02/aa3851.pdf)
- [SOFA docs](https://www.iausofa.org/2021_0512_C/sofa/sofa_pn_c.pdf)
"""
const DCM_ICRF2J2000_BIAS = orient_precession_bias(iau2006a, 0.0)

"""
    DCM_J20002ECLIPJ2000

DCM for the rotation from the Mean Dynamical Equator of J2000 (`MEME2000`) to the 
Mean Ecliptic Equinox. This corresponds to the transformation `J2000 -> ECLIPJ2000` 
in the SPICE toolkit.
"""
const DCM_J20002ECLIPJ2000 = angle_to_dcm(
    orient_obliquity(iau2006a, 0.0), :X
)

"""
    DCM_ICRF2ECLIPJ2000

DCM for the rotation from the International Celestial Reference Frame (`ICRF`) to the 
Mean Ecliptic Equinox of J2000 (`ECLIPJ2000`).
"""
const DCM_ICRF2ECLIPJ2000 = DCM_ICRF2J2000_BIAS * DCM_J20002ECLIPJ2000


# --------------------------------------------------------
# TRANSFORMATIONS
# --------------------------------------------------------

# ICRF -> MEME2000 

function orient_icrf_to_mememod(t::Number)
     # convert TT days since J2000 to TT centuries since J2000
    sec = t/Tempo.CENTURY2DAY

    # fw_angles holds independent on the IAU Model! 
    γ, ϕ, ψ, ε = fw_angles(iau2006b, sec)
    R = fw_matrix(γ, ϕ, ψ, ε)
    return R
end
