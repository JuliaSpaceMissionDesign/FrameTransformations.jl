# [Orient's Public Documentation](@id basic_orient_api)

## [IAU Models](@id iers_models)

This is a list of the supported IAU models and their approximations that can be used to select the desired procedure for the computation of the ITRF to GCRF rotation.

```@docs 

iau2000a
iau2000b

iau2006a
iau2006b 

CPNc
CPNd
```

## IERS Transformations

```@docs 
IERS_EOP
init_eop
prepare_eop
orient_rot3_itrf_to_gcrf
orient_rot6_itrf_to_gcrf
orient_rot9_itrf_to_gcrf
orient_rot12_itrf_to_gcrf
orient_bias_precession
orient_bias_precession_nutation
orient_nutation
orient_obliquity
```

## Moon Transformations

```@docs
orient_rot3_icrf_to_pa421
orient_rot3_icrf_to_pa440
```

## Geodesy 

```@docs 
geoc2pos
pos2geoc
geod2pos
pos2geod
```

## [Axes ID](@id orient_axesid)

This is a list of NAIF IDs for standard axes that are used in astrodynamic applications. 

!!! note 
    Although they are listed in the public documentation section, these IDs are not directly exported by the package.

```@docs 
Orient.AXESID_ICRF
Orient.AXESID_GCRF
Orient.AXESID_ECLIPJ2000
Orient.AXESID_MEME2000
Orient.AXESID_ITRF
Orient.AXESID_MOONME_DE421
Orient.AXESID_MOONPA_DE421
Orient.AXESID_MOONPA_DE440
```

## [Default Rotation Matrices](@id orient_dcms)

```@docs 
Orient.DCM_ICRF_TO_J2000_BIAS
Orient.DCM_ICRF_TO_ECLIPJ2000
Orient.DCM_J2000_TO_ECLIPJ2000
Orient.DCM_MOON_PA421_TO_ME421
Orient.DCM_MOON_PA430_TO_ME421
Orient.DCM_MOON_PA430_TO_ME430
Orient.DCM_MOON_PA440_TO_ME421
```