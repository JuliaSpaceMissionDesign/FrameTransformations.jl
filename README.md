
# FrameTransformations.jl

_A modern high-performance set of tools for transformations between standard and user-defined reference frame._

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaspacemissiondesign.github.io/FrameTransformations.jl/stable/) 
[![Dev Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaspacemissiondesign.github.io/FrameTransformations.jl/dev/) 
[![Build Status](https://github.com/JuliaSpaceMissionDesign/FrameTransformations.jl/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/JuliaSpaceMissionDesign/FrameTransformations.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JuliaSpaceMissionDesign/FrameTransformations.jl/branch/main/graph/badge.svg?token=7fj9BjJhKF)](https://codecov.io/gh/JuliaSpaceMissionDesign/FrameTransformations.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

Are you in search of fundamental routines for efficient and extensible frames transformations?  
If so, this package is the ideal starting point. FrameTransformations.jl is designed to 
provide users with  the ability to create a customized, efficient, flexible, and 
extensible axes/point graph models for mission analysis and space mission design purposes. 

## Features 

- Convert between different time scales and representations (via [`Tempo.jl`](https://github.com/JuliaSpaceMissionDesign/Tempo.jl));
- Read binary ephemeris files (via [`Ephemerides.jl`](https://github.com/JuliaSpaceMissionDesign/Ephemerides.jl));
- Create a custom reference frame systems with user-defined points and axes.
- Transform states between different frames.  

All of this seamlessly integrated with `ForwardDiff.jl`.

## Mission

The development of this package has been performed with the following design goals in mind:

1. **Efficiency**: being a base package a particular attention has been 
    given to the execution time of the different routines as well as most/all of
    them have been optimised and deeply benchmarked.

2. **Extensibility**: attention has been given also to the definition of the 
    interfaces, which have been kept the most essential possible, in such a way 
    their extension can be performed very easily (also thanks to Julia language itself).

3. **Single Responsability Principle**: The different modules in this package 
    have been organized in such a way they are responsible of bringing only *one* 
    of the desired features. This results in the possibility to extend and maybe, 
    in future, detatch some modules to a different package.

4. **Automatic Differentiation**: seamless integration with `ForwardDiff.jl` is targetted 
    to fully exploit its power in higher-level packages constructed on top of `Basic`.
