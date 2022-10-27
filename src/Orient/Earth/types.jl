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