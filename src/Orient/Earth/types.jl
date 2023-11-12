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
