# Welcome to FrameTransformations.jl!

_A modern, high-performance and comprehensive set of tools for transformations between any standard and user-defined reference frame._

Are you in search of fundamental routines for efficient and extensible frames transformations?  
If so, this package is the ideal starting point. FrameTransformations.jl is designed to 
provide users with  the ability to create a customized, efficient, flexible, and 
extensible axes/point graph models for mission analysis and space mission design purposes. 

## Features 

- Convert between different time scales and representations (via [Tempo.jl](https://github.com/JuliaSpaceMissionDesign/Tempo.jl)).
- Read binary ephemeris files (via [Ephemerides.jl](https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl) or [CalcephEphemeris.jl](https://github.com/JuliaSpaceMissionDesign/CalcephEphemeris.jl) extensions).
- Create custom reference frame systems with both standard and user-defined points, axes and directions.
- Transform states and their higher-order derivatives between different frames (up to jerk).
- Compile transformations into zero-overhead, AD-transparent callables for hot loops via [`compile_rotation`](@ref), [`compile_translation`](@ref) and [`compile_direction`](@ref).

Automatic differentiation is supported through [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl), with reverse-mode AD support via [ChainRulesCore.jl](https://github.com/JuliaDiff/ChainRulesCore.jl) and [Mooncake.jl](https://github.com/compintell/Mooncake.jl) package extensions.

## Installation 

This package can be installed using Julia's package manager: 
```julia
julia> import Pkg

julia> Pkg.add("FrameTransformations.jl");
```