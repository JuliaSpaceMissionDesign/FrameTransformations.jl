# Performance Benchmarks 

A comparative analysis of the execution times of various routines of this package against 
other state-of-the-art packages and libraries is here reported. 

## IERS Rotation Models

The Orient sub-module provides the capability to compute high-accuracy Earth rotation models based on the IERS 2010 conventions. However, differently from traditional packages, it leverages Julia's metaprogramming capabilities to generate optimised functions for the computation of the
lengthy 2000A/B nutation series. A comparison against both [ERFA.jl](https://github.com/JuliaAstro/ERFA.jl) and [SatelliteToolbox.jl](https://github.com/JuliaSpace/SatelliteToolbox.jl) is here reported.

```@raw html 
<p align="center">
<img src="https://github.com/JuliaSpaceMissionDesign/FrameTransformations.jl/assets/85893254/1bbd17b6-c487-487a-8f24-99bc671eeabf" width = "512" /> 
</p>
```

!!! details 
    These time benchmarks have been obtained on an Intel Core i7-6700 @ 3.4 GHz with 16 GB of RAM


## Relative Orientation 

Frames' performance have been tested against both [GODOT](https://godot.io.esa.int/docs) and [SPICE.jl](https://github.com/JuliaAstro/SPICE.jl), two of the most-popular packages that can be used to define a system of reference axes and points and to compute the transformations between them. In particular, a comparative analysis of the execution times for computing the relative orientation between different reference axes has been carried out across four distinct scenarios: 

- **Case 1**: From ICRF to the ECLIPJ2000 frame, a time-fixed rotation. 
- **Case 2**: From the IAU-MOON to the ICRF, a time-dependent rotation. 
- **Case 3**: From the IAU-MOON to the ICRF and subsequently to the IAU_EARTH. 
- **Case 4**: From the IAU-MOON, to the ICRF, to the ECLIPJ2000 and finally a shift to the Sun-Earth rotating two-vector frame. Notice that the Sun and Earth positions are retrieved from binary ephemeris kernels.

```@raw html 
<p align="center">
<img src="https://github.com/JuliaSpaceMissionDesign/FrameTransformations.jl/assets/85893254/c8688165-728b-4053-9912-f5ac973e892d" width = "512" /> 
</p>
```

!!! details
    These time benchmarks have been obtained on an Intel Core i7-12700 @ 4.7 GHz with 32 GB of RAM

This figure also underlines the capability of FrameTransformations (here referred to as `Multiverse`) to use different ephemeris readers as backends within the computational graph.

## Relative States 

!!! note 
    Work in progress