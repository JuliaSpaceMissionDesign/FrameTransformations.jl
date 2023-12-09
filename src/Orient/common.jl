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
    orient_rot3_itrf_to_pef(tt::Number)

Compute the rotation matrix from the International Terrestrial Reference Frame (ITRF) to 
the Pseudo-Earth Fixed Frame at time `tt`, expressed in TT seconds since `J2000`.

This is using IAU-76/FK5 Reduction. This is a polar motion only rotation.
Eq. 3-78, Sec. 3.7.3 of Vallado (2013).

### See also 
See also [`IERS_EOP`](@ref).
"""
function orient_rot3_itrf_to_pef(tt::Number)

    if !IERS_EOP.init 
        throw(
            ErrorException(
                "EOP not initialized. Please run 'init_eop' before using this function."
            )
        )
    end

    ttd = tt / Tempo.DAY2SEC
    xₚ = arcsec2rad(interpolate(IERS_EOP.x_TT, ttd))
    yₚ = arcsec2rad(interpolate(IERS_EOP.y_TT, ttd))

    return DCM(
        1.0,    0.0,   -xₚ,
        0.0,    1.0,    yₚ,
         xₚ,    -yₚ,   1.0
    )'

end

"""
    orient_rot3_pef_to_tod(tt::Number; [m]::IAUModel=iau2006a)

Compute the rotation matrix from the Pseudo-Earth Fixed (PEF) to the True Equator of Date 
(TOD) Frame at time `tt`, expressed in TT seconds since `J2000`.

This is using IAU-76/FK5 Reduction. This is a sidereal time only rotation.
Eq. 3-80, Sec. 3.7.3 of Vallado (2013).

### See also 
See also [`orient_gast`](@ref).
"""
function orient_rot3_pef_to_tod(tt::Number; m::IAUModel=iau2006a)
    # TODO: this should use 1980 convention which is not available!
    GAST = orient_gast(m, tt/Tempo.CENTURY2SEC)
    return angle_to_dcm(-GAST, :Z)
end

"""
    orient_rot6_pef_to_tod(tt::Number; [m]::IAUModel=iau2006a)

Compute the rotation matrix and the derivative of the transformation from the Pseudo-Earth 
Fixed (PEF) to the True Equator of Date (TOD) Frame at time `tt`, expressed in TT seconds 
since `J2000`.

This is using IAU-76/FK5 Reduction. This is a sidereal time only rotation.
Eq. 3-80, Sec. 3.7.3 of Vallado (2013).

### See also 
See also [`orient_rot3_pef_to_tod`](@ref), [`orient_gast`](@ref) and [`IERS_EOP`](@ref).
"""
function orient_rot6_pef_to_tod(tt::Number; m::IAUModel=iau2006a)

    if !IERS_EOP.init 
        throw(
            ErrorException(
                "EOP not initialized. Please run 'init_eop' before using this function."
            )
        )
    end

    R = orient_rot3_pef_to_tod(tt; m=m)

    ttd = tt / Tempo.DAY2SEC
    LOD = 1e-3 * interpolate(IERS_EOP.LOD_TT, ttd)
    Rω = skew( SVector{3}(0.0, 0.0, earth_rotation_rate(LOD)) )

    return R, R*Rω
end

"""
    orient_rot3_tod_to_mod(tt::Number; [m]::IAU2006Model=iau2006a)

Compute the rotation matrix from the True Equator of Date (TOD) to the Mean Equator of Date
(MOD)  Frame at time `tt`, expressed in TT seconds since `J2000`.

This is using IAU-76/FK5 Reduction. This is a nutation only rotation.
Eq. 3-86, Sec. 3.7.3 of Vallado (2013).
"""
function orient_rot3_tod_to_mod(tt::Number; m::IAU2006Model=iau2006a)
    # TODO: this should use 1980 convention which is not available!
    t = tt/Tempo.CENTURY2SEC
    Δψ, Δϵ = orient_nutation(m, t)
    ϵ = orient_obliquity(m, t)
    return angle_to_dcm(-ϵ, Δψ, ϵ+Δϵ, :XZX)
end

# --------------------------------------------------------
# I(G)CRF-based transformations
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
    orient_rot3_icrf_to_tod(tt::Number; [m]::IAU2006Model=iau2006a)

Compute the rotation matrix from the International Celestial Reference Frame (ICRF) to 
the True Equator of Date at time `tt`, expressed in TT seconds since `J2000`.

True Equator of Date is obtained applying frame bias, precession and nutation to the ICRF 
pole and origin.
"""
function orient_rot3_icrf_to_tod(tt::Number; m::IAU2006Model=iau2006a)
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

"""
    orient_rot3_icrf_to_pef(tt::Number; [m]::IAU2006Model=iau2006a)

Compute the rotation matrix from the International Celestial Reference Frame (ICRF) to 
the Pseudo Earth Fixed (PEF) Frame at time `tt`, expressed in TT seconds since `J2000`.

### See also 
See also [`orient_rot3_icrf_to_tod`](@ref), [`orient_rot3_pef_to_tod`](@ref) and 
[`orient_rot6_icrf_to_pef`](@ref).
"""
function orient_rot3_icrf_to_pef(tt::Number; m::IAU2006Model=iau2006a)
    R_icrf2tod = orient_rot3_icrf_to_tod(tt; m=m)
    R_tod2pef = orient_rot3_pef_to_tod(tt; m=m)'
    return R_tod2pef * R_icrf2tod
end

"""
    orient_rot6_icrf_to_pef(tt::Number; [m]::IAU2006Model=iau2006a)

Compute the rotation matrix the derivative of the transformation from the International 
Celestial Reference Frame (ICRF) to the Pseudo Earth Fixed (PEF) Frame at time `tt`, 
expressed in TT seconds since `J2000`.

### See also 
See also [`orient_rot3_icrf_to_pef`](@ref).
"""
function orient_rot6_icrf_to_pef(tt::Number; m::IAU2006Model=iau2006a)
    R_icrf2tod = orient_rot3_icrf_to_tod(tt; m=m)
    R_p2t, Ṙ_p2t = orient_rot6_pef_to_tod(tt; m=m)
    return R_p2t' * R_icrf2tod, Ṙ_p2t' * R_icrf2tod
end

"""
    orient_rot3_mod_to_teme(tt::Number; [m]::IAU2006Model=iau2006a)

Compute the rotation matrix from the Mean Equator of Date (MOD) frame to the True Equator, 
Mean Equinox of date at time `tt`, expressed in TT seconds since `J2000`.

This is implemented with a small angle approx of Eq. 4 of Vallado, "Coordinate Frames of the 
US Space Object Catalogs." 
"""
function orient_rot3_mod_to_teme(tt::Number; m::IAU2006Model=iau2006a)
    t = tt / Tempo.CENTURY2SEC 
    Δψ, Δϵ = orient_nutation(m, t)
    ϵ = orient_obliquity(m, t)

    sδϵ = sin(ϵ + Δϵ)

    return DCM(
        1.,      0.,   -Δψ*sδϵ,
        0.,      1.,       -Δϵ,
        Δψ*sδϵ,  Δϵ,        1.
    )
end

"""
    orient_rot3_icrf_to_tod(tt::Number; [m]::IAUModel=iau2006a)

Compute the rotation matrix from the International Celestial Reference Frame (ICRF) to 
the True Equator, Mean Equinox of date at time `tt`, expressed in TT seconds since `J2000`.

### See also 
See also [`@orient_rot3_mod_to_teme`](@ref) and [`orient_rot3_icrf_to_mod`](@ref).
"""
function orient_rot3_icrf_to_teme(tt::Number;  m::IAU2006Model=iau2006a)
    Ricrf2mod = orient_rot3_icrf_to_mod(tt)
    Rmod2teme = orient_rot3_mod_to_teme(tt; m=m)
    return Ricrf2mod * Rmod2teme
end