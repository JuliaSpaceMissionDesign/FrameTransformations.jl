
using StaticArrays: SA 

# ORIENT Precompilation routines 
PrecompileTools.@setup_workload begin
    i2000models = [iau2000a, iau2000b]
    i2006models = [iau2006a, iau2006b]

    PrecompileTools.@compile_workload begin

        # Precompile Earth/IERS routines
        for iaumod in i2006models
            Orient.orient_obliquity(iaumod, 0.0)
            Orient.orient_bias_precession(iaumod, 0.0)
            Orient.orient_bias_precession_nutation(iaumod, 0.0)

            Orient.orient_rot3_itrf_to_gcrf(iaumod, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            Orient.orient_rot6_itrf_to_gcrf(iaumod, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            Orient.orient_rot9_itrf_to_gcrf(iaumod, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            Orient.orient_rot12_itrf_to_gcrf(iaumod, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        end

        for iaumod in i2000models
            Orient.orient_bias_precession_nutation(iaumod, 0.0)

            Orient.orient_rot3_itrf_to_gcrf(iaumod, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            Orient.orient_rot6_itrf_to_gcrf(iaumod, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            Orient.orient_rot9_itrf_to_gcrf(iaumod, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            Orient.orient_rot12_itrf_to_gcrf(iaumod, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        end

        # Ecliptic routines
        Orient.orient_rot3_icrf_to_mememod(0.0)
    end
end

# FRAMES Precompilation routines 
PrecompileTools.@setup_workload begin
    x12 = rand(12)

    x3s = SA[rand(3)...]
    x6s = SA[rand(6)...]
    x9s = SA[rand(9)...]
    x12s = SA[rand(12)...]

    PrecompileTools.@compile_workload begin

        # Precompile twovectors routines 
        Frames.twovectors_to_dcm(x3s, x3s, :XZ)
        Frames.twovectors_to_δdcm(x6s, x6s, :XZ)
        Frames.twovectors_to_δ²dcm(x9s, x9s, :XZ)
        Frames.twovectors_to_δ³dcm(x12s, x12s, :XZ)

        Frames.twovectors_to_dcm(x12, x12, :XZ)
        Frames.twovectors_to_δdcm(x12, x12, :XZ)
        Frames.twovectors_to_δ²dcm(x12, x12, :XZ)
        Frames.twovectors_to_δ³dcm(x12, x12, :XZ)

        Frames._two_vectors_to_rot6(x6s, x6s, :XZ)
        Frames._two_vectors_to_rot9(x9s, x9s, :XZ)
        Frames._two_vectors_to_rot12(x12s, x12s, :XZ)

        Frames._two_vectors_to_rot6(x12, x12, :XZ)
        Frames._two_vectors_to_rot9(x12, x12, :XZ)
        Frames._two_vectors_to_rot12(x12, x12, :XZ)
    end
end
