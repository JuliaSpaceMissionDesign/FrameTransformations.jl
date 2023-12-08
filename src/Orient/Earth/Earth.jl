""" 
    AXESID_ITRF 
   
NAIF Axes ID for the International Terrestrial Reference Frame (ITRF)

!!! note 
    This ID is based upon the ID used to refer to the ITRF93 in NAIF's high-accuracy 
    Earth rotation model PCK kernels.
"""
const AXESID_ITRF = 3000


include("types.jl")
include("fundamentals.jl")
include("eop.jl")

include("builders.jl")

include("obliquity.jl")
include("nutation.jl")
include("precession.jl")

include("iers.jl")
include("sidereal.jl")
include("additional.jl")