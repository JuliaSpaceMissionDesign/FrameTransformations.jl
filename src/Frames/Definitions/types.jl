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
    MEMEMOD

Singleton instance of the [`MeanEquatorMeanEquinoxOfDate`](@ref).

### References 
- [SOFA C](https://www.iausofa.org/2021_0512_C/sofa/sofa_pn_c.pdf)
"""
const MEMEMOD = MeanEquatorMeanEquinoxOfDate()

#
#   Planets IAU
#

"""
    BodyCentricInertialTrueOfDateFrame

An abstract type representing a True of Date body-centric inertial frame.
This is the parent type for the concrete implementations of the frames.

### Example
To define the concrete implementation, of a body-centric frame a type and a 
sigleton shall be defined:

```julia
struct EarthBodyCentricInertialFrame <: BodyCentricInertialTrueOfDateFrame
    vid::Val{Int}
    function EarthBodyCentricInertialFrame()
        new(Val(399))
    end
end

const BCI_EARTH = EarthBodyCentricInertialFrame()
```

Then transformations to another frame of the graph shall be defined. 
Typically, these transformations are given which respect to the 
[`InternationalCelestialReferenceFrame`](@ref), therefore:

```julia
function Rotation(
    origin::InternationalCelestialReferenceFrame,
    target::EarthBodyCentricInertialFrame,
    e::Epoch
)
    ...
    Rotation(origin, target, R)
end
function Rotation(
    origin:::EarthBodyCentricInertialFrame, 
    target::InternationalCelestialReferenceFrame,
    e::Epoch
)
    inv(Rotation(target, origin, e))
end
```

Finally, the transformation has to be added to the frame graph.

```julia
connect!(ICRF, BCI_EARTH)
```
"""
abstract type BodyCentricInertialTrueOfDateFrame <: AbstractBodyCentricInertialFrame end


"""
    BodyCentricInertial2000Frame

An abstract type representing a True of Date body-centric inertial frame.
This is the parent type for the concrete implementations of the frames.

### Example
To define the concrete implementation, of a body-centric frame a type and a 
sigleton shall be defined:

```julia
struct EarthBodyCentricInertial2000Frame <: BodyCentricInertial2000Frame
    vid::Val{Int}
    function EarthBodyCentricInertial2000Frame()
        new(Val(399))
    end
end

const BCI2000_EARTH = EarthBodyCentricInertial2000Frame()
```

Then transformations to another frame of the graph shall be defined. 
Typically, these transformations are given which respect to the 
[`InternationalCelestialReferenceFrame`](@ref), therefore:

```julia
function Rotation(
    origin::InternationalCelestialReferenceFrame,
    target::EarthBodyCentricInertial2000Frame,
    e::Epoch
)
    ...
    Rotation(origin, target, R)
end
function Rotation(
    origin:::EarthBodyCentricInertial2000Frame, 
    target::InternationalCelestialReferenceFrame,
    e::Epoch
)
    inv(Rotation(target, origin, e))
end
```

Finally, the transformation has to be added to the frame graph.

```julia
connect!(ICRF, BCI2000_EARTH)
```
"""
abstract type BodyCentricInertial2000Frame <: AbstractBodyCentricInertialFrame end

# ------------------------------------------------------------------------------
#                                 ROTATING
# ------------------------------------------------------------------------------