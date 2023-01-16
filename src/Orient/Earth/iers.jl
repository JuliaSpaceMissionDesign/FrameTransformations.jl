export orient_itrf_to_gcrf

# Nominal Earth angular velocity  
const ωₑ = 7.292_115_146_706_979e-5

"""
    earth_rotation_rate()
    
Compute the nominal Earth angular velocity. 
"""
earth_rotation_rate() = earth_rotation_rate(0.0)

""" 
    earth_rotation_rate(LOD::Number)

Compute the true angular velocity of the Earth accounting for the 
Length of the Day, i.e., the instantaneous rate of change of UT1 
with respect to a uniform time scale. 
"""
function earth_rotation_rate(LOD::Number)
    return ωₑ*(1-LOD/86400)
end

""" 
    polar_motion(xₚ::N, yₚ::N, t::N, sp::Number)

Compute the Polar Motion rotation matrix from ITRF to TIRS, according to the 
IERS 2010 Conventions.

### Inputs

- `xp`, `yp` -- Coordinates in radians of the Celestial Intermediate Pole 
                (CIP), with respect to the International Terrestrial Reference
                Frame (ITRF).
                
- `t`  -- Terrestrial Time `TT` in Julian centuries since J2000

- `sp` -- The Terrestrial Intermediate Origin (TIO) locator, in radians. It provides 
         the position of the TIO on the equatior fo the CIP. 

### References 
- Luzum, B. and Petit G. (2012), _The IERS Conventions (2010)_, 
[IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
"""
function polar_motion(xₚ::Number, yₚ::Number, t::Number)
    sp = tio_locator(t)
    polar_motion(xₚ, yₚ, t, sp)
end

function polar_motion(xₚ::Number, yₚ::Number, ::Number, sp::Number)
    angle_to_dcm(yₚ, xₚ, -sp, :XYZ)
end


"""
    tio_locator(t::N)

Compute the TIO locator `s'` at date, positioning the Terrestrial Intermediate Origin on 
the equator of the Celestial Intermediate Pole (CIP).

### Input
- `t` -- Terrestrial Time `TT` in Julian centuries since J2000

### Notes 
This function approximates the unpredictable motion of the TIO locator s' with 
its secular drift of ~0.47 μas/century. 

### References
- Luzum, B. and Petit G. (2012), _The IERS Conventions (2010)_, 
[IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- Lambert, S. and Bizouard C. (2002), _Positioning the Terrestrial Ephemeris Origin
in the Terrestrial Reference Frame_, [DOI: 10.1051/0004-6361:20021139](https://www.aanda.org/articles/aa/pdf/2002/40/aa2747.pdf)
"""
function tio_locator(t::Number)
    return -47e-6*t |> arcsec2rad; # arcseconds
end


""" 
    era_rotm(Tᵤ::N)

Compute the TIRS to CIRS Earth Rotation matrix, according to the IERS 2010 
conventions.

### Input 
- `t` -- Julian UT1 date

### References
- Luzum, B. and Petit G. (2012), _The IERS Conventions (2010)_, 
[IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
"""
function era_rotm(t::Number) 
    angle_to_dcm(-earth_rotation_angle(t), :Z)
end


"""
    earth_rotation_angle(t::N) where {N <: Number}

Compute the Earth Rotation Angle (ERA), i.e., the angle between the Celestial 
Intermediate Origin (CIO) and the Terrestrial Intermediate Origin (TIO). 

### Input 
- `t` -- Julian UT1 date

### Notes 
- The function uses the fractional UT1 date to gain additional precision in the 
computations (0.002737.. instead of 1.002737..)

### References 
- Luzum, B. and Petit G. (2012), _The IERS Conventions (2010)_, 
[IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/pfw06.c) library
"""
@inline function earth_rotation_angle(t::Number)
    f = t % 1.0 
    Tᵤ = t - 2451545.0
    return mod2pi(2π * (f + 0.7790572732640 + 0.00273781191135448Tᵤ))
end


"""
    fw2xy(ϵ::Number, ψ::Number, γ::Number, φ::Number)

Compute the CIP X and Y coordinates from Fukushima-Williams bias-precession-nutation 
angles, in radians.

### Inputs 
- `ϵ` -- F-W angle with IAU 2006A/B nutation corrections. 
- `ψ` -- F-W angle with IAU 2006A/B nutation corrections.
- `γ` -- F-W angle  
- `ϕ` -- F-W angle

### References
- Wallace P. T. and Capitaine N. (2006), _Precession-nutation procedures consistent with 
IAU 2006 resolutions_, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/pfw06.c) library
"""
function fw2xy(γ::Number, φ::Number, ψ::Number, ϵ::Number)
    sϵ, cϵ = sincos(ϵ)
    sψ, cψ = sincos(ψ)
    sγ, cγ = sincos(γ)
    sϕ, cϕ = sincos(φ) 

    a = (sϵ*cψ*cϕ - cϵ*sϕ)

    X = sϵ*sψ*cγ - a*sγ
    Y = sϵ*sψ*sγ + a*cγ

    return X, Y
end


"""
    cip_coords(m::IAU2006Model, t::Number)

Computes the CIP X, Y coordinates according to the IAU 2006/2000 A/B. 

### Inputs 
- `m` -- IAU 2006 model 
- `t` -- Terrestrial Time `TT` in Julian Centuries since J2000.0

### References
- Wallace P. T. and Capitaine N. (2006), _Precession-nutation procedures consistent with 
IAU 2006 resolutions_, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/pfw06.c) library
"""
function cip_coords(m::IAU2006Model, t::Number)
 
    # Computes Fukushima-Williams angles
    γ, ϕ, ψ, ϵ = fw_angles(m, t)

    # Computes IAU 2000 nutation components 
    Δψ, Δϵ = orient_nutation(m, t) 
    
    # Applies IAU-2006 compatible nutations 
    ψₙ = ψ + Δψ
    ϵₙ = ϵ + Δϵ
    
    # Retrieves CIP coordinates 
    fw2xy(γ, ϕ, ψₙ, ϵₙ)
end 


include("constants/cio_locator.jl")

# cio_locator(::IAU2006Model, ::Number, ::FundamentalArguments) = ()
build_cio_series(:cio_locator, :IAU2006Model, COEFFS_CIO2006_SP, COEFFS_CIO2006_S)

"""
    cio_locator(m::IAU2006Model, t::Number, x::Number, y::Number)

Compute the CIO Locator `s` in radians, according to the IAU 2010 Conventions.   

### Inputs 
- `m` -- IAU 2006 Model 
- `t` -- Terrestrial Time `TT` in Julian Centuries since J2000.0
- `x, y` -- CIP X, Y coordinates at date `t` in radians.

### References 
- Luzum, B. and Petit G. (2012), _The IERS Conventions (2010)_, 
[IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- Wallace P. T. and Capitaine N. (2006), _Precession-nutation procedures consistent with 
IAU 2006 resolutions_, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/pfw06.c) library
"""
function cio_locator(m::IAU2006Model, t::Number, x::Number, y::Number)
    fa = FundamentalArguments(t, iau2006a)
    cio_locator(m, t, fa)/1e6*π/648000 - x*y/2
end


"""
    cip_motion(m::IAU2006Model, t::Number, dx::Number=0.0, dy::Number=0.0)

Compute the CIRS to GCRS rotation matrix, according to the IAU 2010 Conventions. 

### Inputs 
- `m` -- IAU 2006 Model 
- `t` -- Terrestrial Time `TT` in Julian Centuries since J2000.0
- `dx, dy` -- Optional IERS corrections to account for the free-core nutation and 
time dependent effects. Default are set to zero.

### References 
- Luzum, B. and Petit G. (2012), _The IERS Conventions (2010)_, 
[IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- Wallace P. T. and Capitaine N. (2006), _Precession-nutation procedures consistent with 
IAU 2006 resolutions_, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/pfw06.c) library
"""
function cip_motion(m::IAU2006Model, t::Number, dx::Number=0.0, dy::Number=0.0)

    # Computes CIP X, Y coordinates 
    x, y = cip_coords(m, t)

    x += dx 
    y += dy

    # Computes cio locator `s`
    s = cio_locator(m, t, x, y)
    
    # Retrieves spherical angles E and d. 
    r2 = x^2 + y^2
    E = (r2 > 0) ? atan(y, x) : 0.0 
    d = atan(sqrt(r2 / (1.0 - r2)))

    # This formulation (compared to simplified versions) ensures 
    # that the resulting matrix is orthonormal.
    angle_to_dcm(E+s, -d, -E, :ZYZ)
    
end


"""
    orient_itrf_to_gcrf(t, m::IAU2006Model=iau2006b)

Compute the rotation matrix from ITRF to GCRF at time `t`, according to the IAU 2010 
conventions, using `IAU2006B` as the default model.
"""
function orient_itrf_to_gcrf(t::Number, m::IAU2006Model=iau2006b)
    # TODO: specify time frame!
    # TODO: add computation of xp, yp, dx e dy 
    xp, yp, dx, dy = 0, 0, 0, 0
    return orient_itrf_to_gcrf(t, m)

end

function orient_d_itrf_to_gcrf(t::Number, m::IAU2006Model=iau2006b)

    # TODO: add computation of xp, yp, dx e dy 
    xp, yp, dx, dy = 0, 0, 0, 0
    return orient_d_itrf_to_gcrf(t, m)

end

function orient_dd_itrf_to_gcrf(t::Number, m::IAU2006Model=iau2006b)

    # TODO: add computation of xp, yp, dx e dy 
    xp, yp, dx, dy = 0, 0, 0, 0
    return orient_dd_itrf_to_gcrf(t, m)

end

"""
    orient_itrf_to_gcrf(t, m::IAU2006Model, xₚ, yₚ, δx=0.0, δy=0.0)

Compute the rotation matrix from ITRF to GCRF at time `t`, according to the IAU 2010 
conventions. 
"""
function orient_itrf_to_gcrf(t::Number, m::IAU2006Model, xₚ::Number, yₚ::Number, 
            δx::Number=0.0, δy::Number=0.0)

    # Convert TT since J2000 from days to centuries
    tt_cent = t/Tempo.CENTURY2DAY # FIXME: from days to seconds

    t_ut1 = tt_cent # FIXME: cambia la funzione con quella sotto (una volta creata)
    # t_ut1 = Tempo.tt_to_ut1(t) 

    W = polar_motion(xₚ, yₚ, tt_cent)
    R = era_rotm(t_ut1)
    Q = cip_motion(m, tt_cent, δx, δy)

    return Q*R*W 

end

function orient_d_itrf_to_gcrf(t::Number, m::IAU2006Model, xₚ::Number, yₚ::Number, 
            δx::Number=0.0, δy::Number=0.0)

    # Convert TT since J2000 from days to centuries
    tt_cent = t/Tempo.CENTURY2DAY # FIXME: from days to seconds

    t_ut1 = tt_cent # FIXME: cambia la funzione con quella sotto (una volta creata)
    # t_ut1 = Tempo.tt_to_ut1(t) 

    W = polar_motion(xₚ, yₚ, tt_cent)
    R = era_rotm(t_ut1)
    Q = cip_motion(m, tt_cent, δx, δy)

    ωe = SVector(0.0, 0.0, earth_rotation_rate()) # FIXME: sistema con quella del LOD
    Ω = skew(ωe)

    QR = Q*R

    D = QR*W
    δD = QR*Ω*W

    return D, δD
end


function orient_dd_itrf_to_gcrf(t::Number, m::IAU2006Model, xₚ::Number, yₚ::Number, 
            δx::Number=0.0, δy::Number=0.0)

    # Convert TT since J2000 from days to centuries
    tt_cent = t/Tempo.CENTURY2DAY # FIXME: from days to seconds

    t_ut1 = tt_cent # FIXME: cambia la funzione con quella sotto (una volta creata)
    # t_ut1 = Tempo.tt_to_ut1(t) 

    W = polar_motion(xₚ, yₚ, tt_cent)
    R = era_rotm(t_ut1)
    Q = cip_motion(m, tt_cent, δx, δy)

    ωe = SVector(0.0, 0.0, earth_rotation_rate()) # FIXME: sistema con quella del LOD
    Ω = skew(ωe)

    QR = Q*R 
    QRΩ = QR*Ω

    D = QR*W
    δD = QRΩ*W
    δ²D = QRΩ*Ω*W

    return D, δD, δ²D
end

