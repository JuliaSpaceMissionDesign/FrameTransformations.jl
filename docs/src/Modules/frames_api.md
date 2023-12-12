# [Frames's Public Documentation](@id basic_frames_api)

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
add_axes_bcrtod!
add_axes_computable!
add_axes_eclipj2000!
add_axes_ephemeris!
add_axes_fixedoffset!
add_axes_gcrf!
add_axes_inertial!
add_axes_icrf!
add_axes_itrf!
add_axes_me421!
add_axes_meme2000!
add_axes_mod!
add_axes_tod!
add_axes_teme!
add_axes_pef!
add_axes_pa421!
add_axes_pa440!
add_axes_projected!
add_axes_rotating!
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
Frames.order
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