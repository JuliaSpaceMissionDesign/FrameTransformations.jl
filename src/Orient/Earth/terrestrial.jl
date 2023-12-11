
"""
    orient_rot3_itrf_to_pef(t::Number)

Compute the rotation matrix from the International Terrestrial Reference Frame (ITRF) to 
the Pseudo-Earth Fixed Frame at time `t`, expressed in TT seconds since `J2000`.

### References 
- Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications, Microcosm Press, 
    Hawthorn, CA, USA.

### See also 
See also [`IERS_EOP`](@ref).
"""
function orient_rot3_itrf_to_pef(t::Number)

    check_eop_init()

    ttd = t / Tempo.DAY2SEC

    # Compute the pole coordinates
    xₚ = arcsec2rad(interpolate(IERS_EOP.x_TT, ttd))
    yₚ = arcsec2rad(interpolate(IERS_EOP.y_TT, ttd))

    return angle_to_dcm(xₚ, yₚ, :YX)    

end


"""
    orient_rot3_icrf_to_mod(t::Number, m::IAUModel=iau2006a)

Compute the rotation matrix from the International Celestial Reference Frame (ICRF) to 
the Mean Equinox Mean Equator of Date at time `tt`, expressed in TT seconds since `J2000`.
"""
@inline function orient_rot3_icrf_to_mod(t::Number, m::IAUModel=iau2006a)
    # Compute the bias-precession matrix
    return orient_bias_precession(m, t / Tempo.CENTURY2SEC)
end


function orient_rot3_icrf_to_mod(t::Number, m::IAU1980Model)

    ψₐ, ωₐ, χₐ = precession_angles(m, t)

end