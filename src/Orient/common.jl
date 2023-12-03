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

# function ecliptic_pole(tt::Number)
#     γ, ϕ, _, _ = fw_angles(iau2006b, tt/Tempo.CENTURY2SEC)

#     sγ, cγ = sincos(γ)
#     sϕ, cϕ = sincos(ϕ)
#     return SA[sϕ * sγ, -sϕ*cγ, cϕ]
# end

# function orient_rot3_icrf_to_tod(m::IAU2000Model, tt::Number) # FIXME: not working
#     # Compute ecliptic pole
#     ep = ecliptic_pole(tt)

#     # Compute CIP vector 
#     x, y = cip_coords(m, tt)
#     cip = SA[x, y, sqrt(1 - x^2 - y^2)]

#     v2 = unitvec(cross(cip, ep))
#     v3 = cross(cip, v2)

#     return DCM(v2[1], v2[2], v2[3], v3[1], v3[2], v3[3], cip[1], cip[2], cip[3])
# end

# function orient_rot3_icrf_to_tod(m::IAU2006Model, tt::Number) # FIXME: not working
#     # Compute ecliptic pole
#     γ, ϕ, ψ, ϵ = fw_angles(m, tt/Tempo.CENTURY2SEC)

#     sγ, cγ = sincos(γ)
#     sϕ, cϕ = sincos(ϕ)
#     ep =  SA[sϕ * sγ, -sϕ*cγ, cϕ]

#     # Compute CIP vector
#     Δψ, Δϵ = orient_nutation(m, tt)
#     x, y = fw2xy(γ, ϕ, ψ + Δψ, ϵ + Δϵ)
#     cip = SA[x, y, sqrt(1 - x^2 - y^2)]

#     v2 = unitvec(cross(cip, ep))
#     v3 = cross(cip, v2)

#     return DCM(v2[1], v2[2], v2[3], v3[1], v3[2], v3[3], cip[1], cip[2], cip[3])
# end

# """
#     orient_rot3_icrf_to_tod(tt::Number; iauModel=iau2006b)

# Compute the rotation matrix from the International Celestial Reference Frame (ICRF) to 
# the True Equator of Date at time `tt`, expressed in TT seconds since `J2000`.

# True Equator Of Date is obtained applying frame bias, precession and nutation to the ICRF 
# pole and origin. Based on CIP model and the IAU2006 precession of the ecliptic by default 
# (`iauModel`).
# """
# function orient_rot3_icrf_to_tod(tt::Number; iauModel=iau2006b)
#     # return orient_rot3_icrf_to_tod(iauModel, tt) # FIXME: not working
#     # return orient_bias_precession_nutation(iauModel, tt) # FIXME: not working
# end

# function orient_rot3_icrf_to_teme(m::IAU2006Model, tt::Number)
#     pnb = orient_bias_precession_nutation(m, tt)
#     ee = equinoxes_equation(m, tt)
#     return angle_to_dcm(ee, :Z) * pnb
# end

# """
#     orient_rot3_icrf_to_teme(tt::Number; iauModel=iau2006b)

# Compute the rotation matrix from the International Celestial Reference Frame (ICRF) to 
# the True Equator Mean Equinox at time `tt`, expressed in TT seconds since `J2000`.

# The True Equator Of Date is obtained applying frame bias, precession and nutation to the ICRF 
# pole and origin. These are the TOD axes. To those axes, to recover the Mean Equinox, a rotation 
# about the Z-axis is performed using the equation of equinoxes.

# This transformation is based on CIP model and the IAU2006 precession of the ecliptic by default.
# """
# function orient_rot3_icrf_to_teme(tt::Number; iauModel=iau2006b)
#     return orient_rot3_icrf_to_teme(iauModel, tt)
# end

# function orient_rot3_icrf_to_pef(m::IAU2006Model, tt::Number)
#     pnb = orient_bias_precession_nutation(m, tt)

#     # Compute Earth rotation angle
#     ut1_d = iers_tt_to_ut1(tt * Tempo.CENTURY2DAY)
#     θ = earth_rotation_angle(ut1_d)

#     # Origin equation
#     x, y = bpn2xy(bpn)     
#     s = cio_locator(m, t, x, y)
#     @inbounds begin 
#         ax = x / (1 + bpn[3, 3])
#         xs = 1 - ax*x 
#         ys = - ax * bpn[3, 2]
#         zs = -x 

#         p = bpn[1, 1]*xs + bpn[1, 2]*ys + bpn[1, 3] * zs 
#         q = bpn[2, 1]*xs + bpn[2, 2]*ys + bpn[2, 3] * zs
#     end
#     eo = ((p != 0) || (q != 0)) ? s - atan(q, p) : s

#     # Compute GAST
#     GAST = mod2i(θ - eo)
#     return angle_to_dcm(GAST, :Z) * pnb
# end

# """
#     orient_rot3_icrf_to_pef(tt::Number; iauModel=iau2006b)

# Compute the position rotation matrix from the International Celestial Reference Frame (ICRF) to 
# the Pseudo-Earth Fixed Frame at time `tt`, expressed in TT seconds since `J2000`.

# This transformation is based on CIP model and the IAU2006 precession of the ecliptic by default.
# """
# function orient_rot3_icrf_to_pef(tt::Number; iauModel=iau2006b)
#     return orient_rot3_icrf_to_pef(iauModel, tt)
# end

# """
#     orient_rot6_icrf_to_pef(tt::Number; iauModel=iau2006b)

# Compute the position and velocity rotation matrix from the International Celestial Reference 
# Frame (ICRF) to the Pseudo-Earth Fixed Frame at time `tt`, expressed in TT seconds since `J2000`.

# This transformation is based on CIP model and the IAU2006 precession of the ecliptic by default.
# """
# function orient_rot6_icrf_to_pef(tt::Number; iauModel=iau2006b)
#     R = orient_rot3_icrf_to_pef(iauModel, tt)
#     dθ_dt = 1.00273781191135448 # IERS TN36, Table 1.1 (rev/day) 
#     ω = 2π*dθ_dt/86400.0
#     return R, skew(SA[0, 0, ω])
# end
