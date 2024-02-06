
export DCM_ICRF_TO_MEC2000, DCM_ICRF_TO_EME2000, DCM_EME2000_TO_MEC2000

"""
    DCM_ICRF_TO_EME2000

DCM for the rotation from the International Celestial Reference Frame (`ICRF`) and the 
Mean Equator and Equinox of J2000.0 (`EME2000`). This corresponds to the `J2000` frame in 
the SPICE toolkit.

!!! note 
    The frame bias is here computed using the IAU 2006 Precession model, similarly to ESA's 
    GODOT. Some other software libraries, such as Orekit, use the frame bias of the IAU 2000 
    precession model. The two definitions differ of about 1 arcsecond.

    Moreover, according to [Hilton](https://www.aanda.org/articles/aa/pdf/2004/02/aa3851.pdf) 
    there are multiple possibilities to define the proper rotation between the ICRS and 
    the EME2000. The transformation implemented here correspond to Eq. 6 using the parameters 
    in Table 3, line 1 (RIERS).

### References
- Hilton, James L., and Catherine Y. Hohenkerk. -- Rotation matrix from the mean 
    dynamical equator and equinox at J2000. 0 to the ICRS. -- Astronomy & Astrophysics 
    513.2 (2004): 765-770. DOI: [10.1051/0004-6361:20031552](https://www.aanda.org/articles/aa/pdf/2004/02/aa3851.pdf)
- [SOFA docs](https://www.iausofa.org/2021_0512_C/sofa/sofa_pn_c.pdf)
"""
const DCM_ICRF_TO_EME2000 = IERSConventions.iers_pb(iers2010a, 0.0)

"""
    DCM_EME2000_TO_MEC2000

DCM for the rotation from the Mean Equator and Equinox of J2000 (`EME2000`) to the 
Mean Ecliptic Equinox. This corresponds to the transformation `J2000 -> MEC2000` 
in the SPICE toolkit, and uses the mean obliquity of the ecliptic from the IAU 1976 theory.
"""
const DCM_EME2000_TO_MEC2000 = angle_to_dcm(iers_obliquity(iers1996, 0.0), :X)

"""
    DCM_ICRF_TO_MEC2000

DCM for the rotation from the International Celestial Reference Frame (`ICRF`) to the 
Mean Ecliptic Equinox of J2000 (`MEC2000`).
"""
const DCM_ICRF_TO_MEC2000 = DCM_ICRF_TO_EME2000 * DCM_EME2000_TO_MEC2000
