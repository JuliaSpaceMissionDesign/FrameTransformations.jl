# ------------------------------------------------------------------------------
#                                 INERTIAL
# ------------------------------------------------------------------------------

"""
    InternationalCelestialReferenceFrame 

A type representing the International Celestial Reference Frame.
"""
struct InternationalCelestialReferenceFrame <: AbstractInertialFrame end 

"""
    MeanEclipticOfDate 

A type representing the Mean Ecliptic of Date.
"""
struct MeanEclipticOfDate <: AbstractInertialFrame end 

"""
    MeanEclipticJ2000

A type representing the Mean Ecliptic of J2000 epoch.
"""
struct MeanEclipticJ2000 <: AbstractFrozenFrame end 

# ------------------------------------------------------------------------------
#                                 ROTATING
# ------------------------------------------------------------------------------