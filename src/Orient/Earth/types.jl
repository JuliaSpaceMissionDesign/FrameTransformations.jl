export iau2006a, iau2006b

abstract type IAUModel end
abstract type IAU2006Model <: IAUModel end


########
# 2006 #
########

struct IAU2006A <: IAU2006Model end

"""
    iau2006a

The singleton instance of type `IAU2006a`, representing the IAU 2006A family of models.
"""
const iau2006a = IAU2006A()

struct IAU2006B <: IAU2006Model end

"""
    iau2006a

The singleton instance of type `IAU2006b`, representing the IAU 2006B family of models.
"""
const iau2006b = IAU2006B()

struct IAUSeries{N <: Number}
    sc::N # Sin coefficient
    cc::N # Cos coefficient
    N::SVector{14, Int}
end


"""
    EOPData{T}
EOP Data for IAU 2000A.

!!! note
    Each field will be an `AbstractInterpolation` indexed by the Julian Day.

### Fields
- `x, y`: Polar motion with respect to the crust [arcsec].
- `UT1_UTC`: Irregularities of the rotation angle [s].
- `LOD`: Length of day offset [ms].
- `dX, dY`: Celestial pole offsets referred to the model IAU2000A [milliarcsec].
"""
struct EOPData{T}
    x::T
    y::T
    UT1_UTC::T
    LOD::T
    dX::T
    dY::T
end
