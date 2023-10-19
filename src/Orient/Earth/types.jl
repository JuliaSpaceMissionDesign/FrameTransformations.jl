export iau2006a, iau2006b, iau2000a, iau2000b, iau1980, CPNc, CPNd

abstract type IAUModel end

# ------------
# --- 1980 ---
# ------------
abstract type IAU1980Model <: IAUModel end

struct IAU1980 <: IAU1980Model end

""" 
    iau1980 

The singleton instance of type `IAU1980`, representing the IAU 1980 family of models.
"""
const iau1980 = IAU1980()

# ------------
# --- 2000 ---
# ------------

abstract type IAU2000Model <: IAUModel end

struct IAU2000A <: IAU2000Model end

""" 
    iau2000a 

The singleton instance of type `IAU2006a`, representing the IAU 2000B family of models.
"""
const iau2000a = IAU2000A()

struct IAU2000B <: IAU2000Model end

""" 
    iau2000b 

The singleton instance of type `IAU2006a`, representing the IAU 2000B family of models.
"""
const iau2000b = IAU2000B()

# ------------
# --- 2006 ---
# ------------

abstract type IAU2006Model <: IAUModel end

struct IAU2006A <: IAU2006Model end

"""
    iau2006a

The singleton instance of type `IAU2006a`, representing the IAU 2006A family of models.
"""
const iau2006a = IAU2006A()

struct IAU2006B <: IAU2006Model end

"""
    iau2006b

The singleton instance of type `IAU2006B`, representing the IAU 2006B family of models.

!!! note 
    This is not an official IERS model.
"""
const iau2006b = IAU2006B()

# Approximated concise methods from Capitaine & Wallace: 
# Concise CIO based precession-nutation formulations)
abstract type CPNModel <: IAU2006Model end

struct CPNC <: CPNModel end

"""
    CPNc

The singleton instance of type `CPNC`, representing the concise CPNc from Capitaine & 
Wallace, Concise CIO based precession-nutation formulations, (2008). This model truncates 
the X, Y series to deliver an accuracy of few mas.
"""
const CPNc = CPNC()

struct CPND <: CPNModel end

"""
    CPNd

The singleton instance of type `CPND`, representing the concise CPNd from Capitaine & 
Wallace, Concise CIO based precession-nutation formulations, (2008). This model truncates 
the X, Y series to deliver an accuracy of few arcseconds.
"""
const CPNd = CPND()

struct IAUSeries{N<:Number}
    sc::N # Sin coefficient
    cc::N # Cos coefficient
    N::SVector{14,Int}
end

"""
    EOPData{T}

EOP Data for IAU 2000A.

!!! note
    Each field will be an `AbstractInterpolation` indexed by the Julian Day.

!!! note
    A set of EOP data parameterized as function of the Terrestrial Time (TT) scale is 
    also provided. These leverage the fact that the difference between `TT` and `UT1` 
    is a continuous and differentiable function. Passing through `UTC` would introduce 
    discontinuities across the introduction of leap-seconds, which would in turn cause 
    numerical issues when the function derivatives are accumulated over time.

### Fields
- `x, y`: Polar motion with respect to the crust [arcsec].
- `UT1_UTC`: Irregularities of the rotation angle [s].
- `LOD`: Length of day offset [ms].
- `dX, dY`: Celestial pole offsets referred to the model IAU2000A [milliarcsec].
- `x_TT, y_TT`: Polar motion parameterized by TT days [arcsec].
- `UT1_TT`: Irregularities of the rotation angle parameterized by TT days [s].
- `LOD_TT`: Length of day offset parameterized by TT days [ms].
- `dX_TT, dY_TT`: Celestial pole offsets parameterized by TT days.
"""
struct EOPData{T}
    x::T
    y::T
    UT1_UTC::T
    LOD::T
    dX::T
    dY::T

    # EOP data parametrized by TT epoch
    # These are parsed automatically in eop.jl to allow direct computation of EOP 
    # without performing many transformations from TT/TDB, which are considered to be equal.
    x_TT::T
    y_TT::T
    UT1_TT::T
    LOD_TT::T
    dX_TT::T
    dY_TT::T
end
