module Ephemeris

using Basic.Tempo

import InterfacesUtils.IO: load
import InterfacesUtils.Interfaces
using InterfacesUtils.Interfaces.Ephemeris
using InterfacesUtils.Interfaces.Errors

import PrecompileTools

# include("abstract.jl")
include("empty.jl")
include("calceph.jl")

# Precompilation routines 
PrecompileTools.@setup_workload begin 

    PrecompileTools.@compile_workload begin 

        ephem_timescale(NullEphemerisProvider())
        
    end
end

end
