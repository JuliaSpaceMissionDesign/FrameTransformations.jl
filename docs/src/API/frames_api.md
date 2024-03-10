# [Public Documentation](@id basic_frames_api)

## Frame System

```@docs 
FrameSystem
```

## Axes 

```@docs
@axes
axes_alias

is_inertial
is_timefixed

ComputableAxesVector

add_axes_bci2000!
add_axes_bcifix!
add_axes_bcrtod!
add_axes_cirf!
add_axes_computable!
add_axes_ecl2000!
add_axes_eme2000!
add_axes_ephemeris!
add_axes_fixedoffset!
add_axes_gcrf!
add_axes_gtod!
add_axes_inertial!
add_axes_icrf!
add_axes_itrf!
add_axes_me421!
add_axes_mod!
add_axes_pef!
add_axes_pa421!
add_axes_pa440!
add_axes_projected!
add_axes_rotating!
add_axes_tirf!
add_axes_tod!
add_axes_topocentric!

```

## Points

```@docs
@point
point_alias

add_point_dynamical!
add_point_ephemeris!
add_point_fixed!
add_point_root!
add_point_surface!
add_point_updatable!

update_point!
```

## [Rotations](@id rotation_api)

```@docs 
Rotation
FrameTransformations.order
Base.inv

rotation3
rotation6
rotation9
rotation12
```

## Transformations 

```@docs 

LightTime
PlanetaryAberration

vector3
vector6
vector9 
vector12 
```

## Axes ID

This is a list of NAIF IDs for standard axes that are used in astrodynamic applications.

!!! note 
    Although they are listed in the public documentation section, these IDs are not directly exported by the package.

```@docs 
FrameTransformations.AXESID_ICRF
FrameTransformations.AXESID_MOONME_DE421
FrameTransformations.AXESID_MOONPA_DE421
FrameTransformations.AXESID_MOONPA_DE440
FrameTransformations.AXESID_ECL2000
FrameTransformations.AXESID_EME2000
FrameTransformations.AXESID_GCRF
```

## Default Rotation Matrices 

```@docs 
DCM_ICRF_TO_EME2000
DCM_MOON_PA421_TO_ME421
DCM_MOON_PA430_TO_ME421
DCM_MOON_PA430_TO_ME430
DCM_MOON_PA440_TO_ME421
```
