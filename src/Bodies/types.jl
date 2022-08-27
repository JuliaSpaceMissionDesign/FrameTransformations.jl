export CelestialBody,
       Barycenter,
       Planet,
       NaturalSatellite, 
       MinorBody,
       Asteroid, 
       Comet,

       NAIFId,

       body_declination,
       body_declination_rate,
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
    body_from_naifid(id::NAIFId)

Return a celestial body instance based on its NAIF ID code.
"""
body_from_naifid(id::NAIFId) = body_from_naifid(Val(id))

"""
    body_parent(body::CelestialBody)::NAIFId
    body_parent(body::NAIFId)::NAIFId

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
    body_gm(body::NAIFId)::Float64

Return the gravitational parameter ``\\mu = GM`` of `body` in km^3/s^2.

# References
- [NASA NAIF](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc)
"""
function body_gm end

"""
    body_mean_radius(body::CelestialBody)::Float64
    body_mean_radius(body::NAIFId)::Float64

Return the mean radius of `body` in km.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)
"""
function body_mean_radius end

"""
    body_polar_radius(body::CelestialBody)::Float64
    body_polar_radius(body::NAIFId)::Float64

Return the polar radius of `body` in km.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)
"""
function body_polar_radius end

"""
    body_equatorial_radius(body::CelestialBody)::Float64
    body_equatorial_radius(body::NAIFId)::Float64

Return the polar radius of `body` in km.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)
"""
function body_equatorial_radius end

for fun in (:body_parent, :body_system_equivalent, :body_equatorial_radius,
    :body_polar_radius, :body_mean_radius, :body_gm)
    @eval begin
        $fun(id::NAIFId) = $fun(Val(id))
    end
end
