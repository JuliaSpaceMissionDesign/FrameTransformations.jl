"""
    ecliptic_pole(m::IAU2006Model, t::Number)

Computes ecliptic pole of date in ICRF.

- Wallace, P. T., & Capitaine, N. (2006). Precession-nutation procedures consistent with IAU 
2006 resolutions. Astronomy & Astrophysics, Eq. 20.
"""
function ecliptic_pole(m::IAU2006Model, t::Number)
    γ, ϕ, _, _ = fw_angles(m, t)
    return SVector{3}(sin(ϕ)*sin(γ), -sin(ϕ)*cos(γ), cos(ϕ))
end

"""
    cip(m::IAUModel, t::Number)

Compute Celestial Intermediate Pole vector.

- Wallace, P. T., & Capitaine, N. (2006). Precession-nutation procedures consistent with IAU 
2006 resolutions. Astronomy & Astrophysics.
"""
function cip(m::IAUModel, t::Number)
    xs, ys = cip_coords(m, t)
    return SVector{3}(xs, ys, sqrt(1 -(xs^2 + ys^2)))
end