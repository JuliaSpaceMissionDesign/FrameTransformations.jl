export ICRF, 
       MEME2000,
       ECLIPJ2000,
       MEMEOD,
       
       # types 
       InternationalCelestialReferenceFrame, 
       MeanEquatorMeanEquinoxJ2000,
       EclipticEquinoxJ2000,
       MeanEquatorMeanEquinoxOfDate

# ------------------------------------------------------------------------------
#                                 INERTIAL
# ------------------------------------------------------------------------------
#
#   Celestial
#

"""
    InternationalCelestialReferenceFrame 

A type representing the International Celestial Reference Frame.
"""
struct InternationalCelestialReferenceFrame <: AbstractInertialFrame end 

"""
    ICRF

Singleton instance of the [`InternationalCelestialReferenceFrame`](@ref).
"""
const ICRF = InternationalCelestialReferenceFrame()
register!(ICRF)

"""
    MeanEquatorMeanEquinoxJ2000

A type representing the Mean Equator Mean Equinox of J2000.
"""
struct MeanEquatorMeanEquinoxJ2000 <: AbstractInertialFrame end

"""
    MEME2000

Singleton instance of the [`MeanEquatorMeanEquinoxJ2000`](@ref). This is the 
realization of the spice `J2000` frame.

### References 
- [SPICE](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/Tutorials/pdf/individual_docs/17_frames_and_coordinate_systems.pdf) documentation.
"""
const MEME2000 = MeanEquatorMeanEquinoxJ2000()

"""
    EclipticEquinoxJ2000

A type representing the Ecliptic and Equinox of J2000.
"""
struct EclipticEquinoxJ2000 <: AbstractInertialFrame end

"""
    ECLIPJ2000

Singleton instance of the [`EclipticEquinoxJ2000`](@ref).
"""
const ECLIPJ2000 = EclipticEquinoxJ2000()

"""
    MeanEquatorMeanEquinoxOfDate

A type representing the Mean Ecliptic and Equinox of Date.
"""
struct MeanEquatorMeanEquinoxOfDate <: AbstractInertialFrame end

"""
    MEMEOD

Singleton instance of the [`MeanEquatorMeanEquinoxOfDate`](@ref).
"""
const MEMEOD = MeanEquatorMeanEquinoxOfDate()

# ------------------------------------------------------------------------------
#                                 ROTATING
# ------------------------------------------------------------------------------