export CelestialBody,
       Barycenter,
       Planet,
       NaturalSatellite, 
       MinorBody,
       Asteroid, 
       Comet,

       NAIFId,

       body_along_orbit_radius,
       body_declination,
       body_declination_rate,
       body_ellipsoid,
       body_equatorial_radius,
       body_from_naifid,
       body_gm,
       body_mean_radius,
       body_naifid,
       body_parent,
       body_polar_radius,
       body_right_ascension,
       body_right_ascension_rate,
       body_rotation_angle,
       body_rotation_rate,
       body_subplanetary_radius,
       body_system_equivalent

"""
NAIFId
   
   An integer code to identify celestial bodies and other objects in space.
   
   # References
   - [NASA NAIF](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html)
"""
const NAIFId = Integer

"""
    AbstractBody

Abstract supertype for all bodies.
"""
abstract type AbstractBody end

"""
    CelestialBody

Abstract supertype for all celestial bodies and pseudo-bodies.
"""
abstract type CelestialBody <: AbstractBody end

"""
    Barycenter <: CelestialBody

Abstract supertype for representing the barycenters of the solar system 
and planetary systems as pseudo-bodies.
"""
abstract type Barycenter <: CelestialBody end

"""
    Planet <: CelestialBody

Abstract supertype for planets.
"""
abstract type Planet <: CelestialBody end

"""
    NaturalSatellite <: CelestialBody

Abstract supertype for all natural satellites (moons).
"""
abstract type NaturalSatellite <: CelestialBody end

"""
    MinorBody <: CelestialBody

Abstract supertype for minor solar system bodies.
"""
abstract type MinorBody <: CelestialBody end

"""
    Asteroid <: MinorBody

Abstract supertype for asteroids.
"""
abstract type Asteroid <: MinorBody end

"""
    Comet <: MinorBody

Abstract supertype for comets.
"""
abstract type Comet <: MinorBody end

"""
    to_naifid(body::CelestialBody)::NAIFId

Get the NAIF ID code for `body`.
"""
function body_naifid end

"""
    from_naifid(id::NAIFId)

Return a celestial body instance based on its NAIF ID code.
"""
body_from_naifid(id::NAIFId) = body_from_naifid(Val(id))

"""
    body_parent(body::CelestialBody)::NAIFId

Get parent of a given body
"""
function body_parent end 

"""
    body_system_equivalent(body::CelestialBody)::NAIFId

Return the body system equivalent body or barycenter.
"""
function body_system_equivalent end 

"""
    body_gm(body::CelestialBody)::Float64

Return the gravitational parameter ``\\mu = GM`` of `body` in km^3/s^2.

# References
- [NASA NAIF](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc)
"""
function body_gm end

"""
    body_mean_radius(body::CelestialBody)::Float64

Return the mean radius of `body` in km.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)
"""
function body_mean_radius end

"""
    body_polar_radius(body::CelestialBody)::Float64

Return the polar radius of `body` in km.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)
"""
function body_polar_radius end

"""
    body_equatorial_radius(body::CelestialBody)::Float64

Return the polar radius of `body` in km.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)
"""
function body_equatorial_radius end

"""
    body_subplanetary_radius(body::CelestialBody)::Float64

Return the polar radius of `body` in km.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)
"""
function body_subplanetary_radius end

"""
    body_along_orbit_radius(body::CelestialBody)::Float64

Return the polar radius of `body` in km.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)
"""
function body_along_orbit_radius end

"""
    ellipsoid(body::CelestialBody)::Float64

Return the subplanetary, along-orbit, and polar radii of `body` in km which constitute its
tri-axial ellipsoid.

# Example
```jldoctest
julia> ellipsoid(earth)
(6378.1366, 6378.1366, 6356.7519)
```
# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)
"""
function body_ellipsoid(body::CelestialBody) where {T}
    return body_subplanetary_radius(body), body_along_orbit_radius(body), body_polar_radius(body)
end

"""
    right_ascension(body::CelestialBody, ep::Float64)

Return the right ascension of the body-fixed frame of `body` w.r.t. the ICRF 
at epoch `ep` in rad. Here `ep` is expressed in centuries since J2000 in the 
TDB time scale.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)
"""
function body_right_ascension end

"""
    right_ascension_rate(body::CelestialBody, ep::Float64)

Return the right ascension rate of change of the body-fixed frame of `body` w.r.t. the
ICRF at epoch `ep` in rad/s. Here `ep` is expressed in centuries since J2000 in the 
TDB time scale.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)
"""
function body_right_ascension_rate end

"""
    body_declination(body::CelestialBody, ep::Float64)

Return the declination of the body-fixed frame of `body` w.r.t. the ICRF 
at epoch `ep` in rad. Here `ep` is expressed in centuries since J2000 in the 
TDB time scale.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)
"""
function body_declination end

"""
    body_declination_rate(body::CelestialBody, ep::Float64)

Return the declination rate of change of the body-fixed frame of `body` w.r.t. the
ICRF at epoch `ep` in rad/s. Here `ep` is expressed in centuries since J2000 in the 
TDB time scale.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)
"""
function body_declination_rate end

"""
    body_rotation_angle(body::CelestialBody, ep::Float64)

Return the declination of the body-fixed frame of `body` w.r.t. the ICRF 
at epoch `ep` in rad. Here `ep` is expressed in centuries since J2000 in the 
TDB time scale.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)
"""
function body_rotation_angle end

"""
    body_rotation_rate(body::CelestialBody, ep::Float64)

Return the declination rate of change of the body-fixed frame of `body` w.r.t. the
ICRF at epoch `ep` in rad/s. Here `ep` is expressed in centuries since J2000 in the 
TDB time scale.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)
"""
function body_rotation_rate end