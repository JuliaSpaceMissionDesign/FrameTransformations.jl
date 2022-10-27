"""
    orient_body_angles(b::Union{NAIFId, B}, t::N) where {N<:AbstractFloat, B <: CelestialBody}

Construct Euler angles for planetary reference frames definition. 

### Inputs 
- `b` -- body name 
- `t` -- `TDB` centuries since J2000

### Outputs 
Euler angles:
- `90 + α` -- node axis location (intersection between the body and the ICRF equator)
- `90 - δ` -- inclination of the planet’s equator to the celestial equator
- `W` -- position of the prime meridian 

### References
- Archinal, B. A., et al. "Report of the IAU working group on cartographic coordinates 
    and rotational elements: 2015." Celestial Mechanics and Dynamical Astronomy 
    130.3 (2018): 1-46, [doi](https://link.springer.com/content/pdf/10.1007/s10569-017-9805-5.pdf)
"""
function orient_body_angles(b::Union{NAIFId, B}, t::N) where {N<:AbstractFloat, B <: CelestialBody}
    orient_right_ascension(b, t) + π/2, π/2 - orient_declination(b, t), mod2pi(orient_rotation_angle(b, t))
end

"""
    orient_body_rates(b::Union{NAIFId, B}, t::N) where {N<:AbstractFloat, B <: CelestialBody}

Construct Euler angles rates for planetary reference frames definition. 

### Inputs 
- `b` -- body name 
- `t` -- `TDB` centuries since J2000

### Outputs 
Rates of the Euler angles defined as:
- `90 + α` -- node axis location (intersection between the body and the ICRF equator)
- `90 - δ` -- inclination of the planet’s equator to the celestial equator
- `W` -- position of the prime meridian 

### References
- Archinal, B. A., et al. "Report of the IAU working group on cartographic coordinates 
    and rotational elements: 2015." Celestial Mechanics and Dynamical Astronomy 
    130.3 (2018): 1-46, [doi](https://link.springer.com/content/pdf/10.1007/s10569-017-9805-5.pdf)
"""
function orient_body_rates(b::Union{NAIFId, B}, t::N) where {N<:AbstractFloat, B <: CelestialBody}
    orient_right_ascension_rate(b, t), -orient_declination_rate(b, t), orient_rotation_rate(b, t)
end