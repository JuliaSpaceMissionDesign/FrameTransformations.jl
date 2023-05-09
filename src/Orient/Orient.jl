module Orient

using DelimitedFiles
using LinearAlgebra
using Logging
using ReferenceFrameRotations
using RemoteFiles: @RemoteFile, download, path
using StaticArrays

using SMDInterfacesUtils.Interfaces.Ephemeris: AbstractEphemerisProvider, 
                                                ephem_available_axes, 
                                                ephem_orient!
                                                
using Tempo
using FrameTransformations.Utils: skew

using SMDInterfacesUtils.Math: InterpAkima, interpolate, arcsec2rad

using PrecompileTools: PrecompileTools

# Earth
include("Earth/Earth.jl")

# Moon 
include("moon.jl")

# Planets
include("planets.jl")

# Ecliptic 
include("ecliptic.jl")

# Precompilation routines 
PrecompileTools.@setup_workload begin
    i2000models = [iau2000a, iau2000b]
    i2006models = [iau2006a, iau2006b]

    PrecompileTools.@compile_workload begin

        # Precompile Earth/IERS routines
        for iaumod in i2006models
            orient_obliquity(iaumod, 0.0)
            orient_bias_precession(iaumod, 0.0)
            orient_bias_precession_nutation(iaumod, 0.0)

            orient_rot3_itrf_to_gcrf(iaumod, 0.0)
            orient_rot6_itrf_to_gcrf(iaumod, 0.0)
            orient_rot9_itrf_to_gcrf(iaumod, 0.0)
            orient_rot12_itrf_to_gcrf(iaumod, 0.0)
        end

        for iaumod in i2000models
            orient_bias_precession_nutation(iaumod, 0.0)

            orient_rot3_itrf_to_gcrf(iaumod, 0.0)
            orient_rot6_itrf_to_gcrf(iaumod, 0.0)
            orient_rot9_itrf_to_gcrf(iaumod, 0.0)
            orient_rot12_itrf_to_gcrf(iaumod, 0.0)
        end

        # Ecliptic routines
        orient_rot3_icrf_to_mememod(0.0)
    end
end

end
