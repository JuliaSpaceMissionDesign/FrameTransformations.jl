export DCM_ICRF_TO_MEME2000

"""
    AXESID_MEME2000

Axes ID for the Mean Dynamical Equator and Equinox of J2000.0. 

!!! note 
    In SPICE the J2000 and ICRF axes are considered equal, thus there exist no 
    specific NAIF ID for the MEME2000 axes. 22 has been chosen because it is the 
    first unassigned axes ID among the built-in SPICE frames. 
"""
const AXESID_MEME2000 = 22

"""
    DCM_ICRF_TO_MEME2000

DCM for the rotation from the International Celestial Reference Frame (`ICRF`) and the 
Mean Equator and Equinox of J2000.0 (`MEME2000`). This corresponds to the `J2000` frame in 
the SPICE toolkit.

!!! note 
    The frame bias is here computed using the IAU 2006 Precession model, similarly to ESA's 
    GODOT. Some other software libraries, such as Orekit, use the frame bias of the IAU 2000 
    precession model. The two definitions differ of about 1 arcsecond.

    Moreover, according to [Hilton](https://www.aanda.org/articles/aa/pdf/2004/02/aa3851.pdf) 
    there are multiple possibilities to define the proper rotation between the ICRS and 
    the MEME2000. The transformation implemented here correspond to Eq. 6 using the parameters 
    in Table 3, line 1 (RIERS).

### References
- Hilton, James L., and Catherine Y. Hohenkerk. -- Rotation matrix from the mean 
    dynamical equator and equinox at J2000. 0 to the ICRS. -- Astronomy & Astrophysics 
    513.2 (2004): 765-770. DOI: [10.1051/0004-6361:20031552](https://www.aanda.org/articles/aa/pdf/2004/02/aa3851.pdf)
- [SOFA docs](https://www.iausofa.org/2021_0512_C/sofa/sofa_pn_c.pdf)
"""
const DCM_ICRF_TO_MEME2000 = orient_bias_precession(iau2006a, 0.0)

# --------------------------------------------------------
# TRANSFORMATIONS
# --------------------------------------------------------

"""
    orient_rot3_icrf_to_mod(tt::Number)

Compute the rotation matrix from the International Celestial Reference Frame (ICRF) to 
the Mean Equinox Mean Equator of Date at time `tt`, expressed in TT seconds since `J2000`.

Mean Equator Of Date is obtained applying frame bias and precession to the ICRF pole and origin.
Fukushima-Williams parametrization for the equator and ecliptic precession is used. 
Consistent with the IAU2006 precession model.
"""
function orient_rot3_icrf_to_mod(tt::Number)
    # convert TT seconds since J2000 to TT centuries since J2000
    T = tt / Tempo.CENTURY2SEC

    # fw_angles holds independent on the IAU Model! 
    γ, ϕ, ψ, ε = fw_angles(iau2006b, T)
    R = fw_matrix(γ, ϕ, ψ, ε)
    return R
end

"""
    orient_rot3_icrf_to_tod(tt::Number; [m]::IAUModel=iau2006a)

Compute the rotation matrix from the International Celestial Reference Frame (ICRF) to 
the True Equator of Date at time `tt`, expressed in TT seconds since `J2000`.

True Equator of Date is obtained applying frame bias, precession and nutation to the ICRF 
pole and origin.
"""
function orient_rot3_icrf_to_tod(tt::Number; m::IAUModel=iau2006a)
    t = tt / Tempo.CENTURY2SEC

    # Compute CIP vector 
    xs, ys = cip_coords(m, t)
    Ĉ = SVector{3}(xs, ys, sqrt(1 -(xs^2 + ys^2)))

    # Compute ecliptic pole 
    K = ecliptic_pole(iau2006a, t)
    
    # Compute rotation matrix 
    X̂ = unitvec(cross(Ĉ, K))
    return DCM(hcat(X̂, cross(Ĉ, X̂), Ĉ)')

end