# [Light Time Corrections](@id tutorial_04_lighttime)

In this tutorial we will explore FrameTransformations.jl capabilities to account for aberration corrections when computing the relative position and velocity between two points. Indeed, this package is capable of replicating all of SPICE's aberration correction features. 

Assuming we are in the **Reception** case, in which photons depart from the target's location at the light-time correted epoch `et-lt` and arrive at the observer's location at `et`, the following corrections are available:
- **One-way Light Time**: it provides the position of the target at the moment it emitted photons arriving at the observer at `et`.
- **Stellar Aberration**: it corrects for one-way light time aswell as stellar aberration by modifiying the relative position to account for the observer's velocity relative to the solar system barycenter. The output is the apparent position of the target as seen by the observer.

For a more detailed overview of this concepts, please refer to the original [SPICE documentation](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkezp_c.html).

## One-way Light Time 

!!! note 
    Work in progress

## Stellar Aberration 


!!! note 
    Work in progress