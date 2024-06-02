# [Axes](@id axes_api) 

## Core

```@docs
add_axes!

add_axes_root!
add_axes_inertial!
add_axes_rotating!
add_axes_fixedoffset!
add_axes_alias!
add_axes_frozen!
add_axes_ephemeris!
```

## Celestial

```@docs
add_axes_gcrf!
add_axes_icrf!
add_axes_eme2000!
```

## Ecliptic

```@docs
add_axes_ecl2000!
```

## Terrestrial 

```@docs
add_axes_itrf!
add_axes_tirf!
add_axes_cirf!

add_axes_mod!
add_axes_tod!
add_axes_gtod!
add_axes_pef!
```

## Planetary

```@docs
add_axes_bci2000!
add_axes_bcrtod!
```

## Lunar

```@docs
add_axes_me421!
add_axes_pa421!
add_axes_pa440!
```

## Topocentric 

```@docs
add_axes_topocentric!
```

## Others 

```@docs
add_axes_twovectors!
```

## Utils

### [IDs](@id frames_axesid)

This is a list of NAIF IDs for standard axes that are used in astrodynamic applications.

```@docs 
FrameTransformations.AXESID_ICRF
FrameTransformations.AXESID_MOONME_DE421
FrameTransformations.AXESID_MOONPA_DE421
FrameTransformations.AXESID_MOONPA_DE440
FrameTransformations.AXESID_ECL2000
FrameTransformations.AXESID_EME2000
FrameTransformations.AXESID_GCRF
```

### [Rotation Matrices](@id frames_dcms)

```@docs 
FrameTransformations.DCM_ICRF_TO_EME2000
```