export CelestialBody,
       Barycenter,
       Planet,
       NaturalSatellite, 
       MinorBody,
       Asteroid, 
       Comet,

       BodyId,

       body_equatorial_radius,
       body_gm,
       body_mean_radius,
       body_naifid,
       body_parent,
       body_polar_radius,
       body_system_equivalent

"""
    BodyId
   
An integer code to identify celestial bodies and other objects in space.
   
# References
- [NASA NAIF](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html)
"""
struct BodyId
    astroid::Int
    naifid::Int
    BodyId(astroid::Integer, naifid::Integer) = new(astroid, naifid)
end

function Base.show(io::IO, id::BodyId)
    println(io, "BodyId(astroid=$(id.astroid), naifid=$(id.naifid))")
end

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
