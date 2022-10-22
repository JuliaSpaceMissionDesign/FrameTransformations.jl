export orient_icrs2ecleod

"""
    orient_icrs2ecleod(jd1tt::N, jd2tt::N2) where {N<:Number, N2<::Number}

Rotation matrix from ICRS equatorial to ecliptic and equinox of date, IAU 2006.

### Input/s 
- `jd1tt`, `jd2tt` -- TT as a 2-part Julian Date

### Output 
Rotation matrix from ICRS to ecliptic.

The matrix is in the sense E_ep = rm x P_ICRS, where P_ICRS is a vector with 
respect to ICRS right ascension and declination axes and E_ep is the same vector 
with respect to the (inertial) ecliptic and equinox of date.

P_ICRS is a free vector, merely a direction, typically of unit magnitude, and 
not bound to any particular spatial origin, such as the Earth, Sun or SSB.  
No assumptions are made about whether it represents starlight and embodies 
astrometric effects such as parallax or aberration.  The transformation is 
approximately that between mean J2000.0 right ascension and declination and ecliptic
longitude and latitude, with only frame bias (always less than 25 mas) to disturb 
this classical picture.
"""
function orient_icrs2ecleod(jd1tt::N, jd2tt::N2) where {N<:Number, N2<:Number}
    # Compute obliquity
    ε = obliquity(iau2006, jd1tt, jd2tt)

    # Compute precession bias matrix
    bp = precession_bias(iau2006, jd1tt, jd2tt)

    # Equatorial of date to ecliptic matrix.
    R_eod_ecl = angle_to_dcm(ε, :X)

    # ICRS to ecliptic coordinates rotation matrix, IAU 2006
    return R_eod_ecl * bp
end