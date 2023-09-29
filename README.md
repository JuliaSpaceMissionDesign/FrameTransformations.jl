
# FrameTransformations.jl

_A modern high-performance set of tools for transformations between standard and user-defined reference frame._

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaspacemissiondesign.github.io/FrameTransformations.jl/stable/) 
[![Dev Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaspacemissiondesign.github.io/FrameTransformations.jl/dev/) 
[![Build Status](https://github.com/JuliaSpaceMissionDesign/FrameTransformations.jl/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/JuliaSpaceMissionDesign/FrameTransformations.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JuliaSpaceMissionDesign/FrameTransformations.jl/branch/main/graph/badge.svg?token=7fj9BjJhKF)](https://codecov.io/gh/JuliaSpaceMissionDesign/FrameTransformations.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

Are you in search of fundamental routines for efficient and extensible frames transformations?  If so, `FrameTransformations` is the ideal starting point. The package is designed to provide users with the ability to create a customized, efficient, flexible, and extensible axes/point graph models for mission analysis and space mission design purposes. 

## Installation 

This package can be installed using Julia's package manager: 
```julia 
julia> import Pkg 

julia> Pkg.add("FrameTransformations.jl");
```

## Features 

- Convert between different time scales and representations (via [Tempo.jl](https://github.com/JuliaSpaceMissionDesign/Tempo.jl));
- Read binary ephemeris files (via [Ephemerides.jl](https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl) or [CalcephEphemeris.jl](https://github.com/JuliaSpaceMissionDesign/CalcephEphemeris.jl))
- Create custom reference frame systems with both standard and user-defined points and axes.
- Transform states and their higher-order derivatives between different frames (up to jerk)

All of this seamlessly integrated with [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl).

## Documentation 
For further information on this package and its tutorials please refer to the 
[stable documentation](https://juliaspacemissiondesign.github.io/FrameTransformations.jl/stable/).

## Support
If you found this package useful, please consider starring the repository. We also encourage 
you to take a look at other astrodynamical packages of the [JSMD](https://github.com/JuliaSpaceMissionDesign/) organisation.