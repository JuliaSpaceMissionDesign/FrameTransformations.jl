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

## EOP Data 

```@docs
init_eop
prepare_eop
eop_data_filename 
IERS_EOP
```

## IERS Transformations

```@docs
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
Orient.AXESID_B1950
Orient.AXESID_ECLIPB1950
Orient.AXESID_ECLIPJ2000
Orient.AXESID_FK4
Orient.AXESID_GALACTIC 
Orient.AXESID_GCRF
Orient.AXESID_ICRF
Orient.AXESID_ITRF
Orient.AXESID_MEME2000
Orient.AXESID_MOONME_DE421
Orient.AXESID_MOONPA_DE421
Orient.AXESID_MOONPA_DE440
```

## [Default Rotation Matrices](@id orient_dcms)

!!! note 
    Although they are listed in the public documentation section, the rotation matrices of some older frames (e.g., B1950, FK4 and GALACTIC) are not exported by the package.

```@docs 
Orient.DCM_B1950_TO_ECLIPB1950
Orient.DCM_B1950_TO_FK4
Orient.DCM_ICRF_TO_ECLIPJ2000
Orient.DCM_ICRF_TO_MEME2000
Orient.DCM_FK4_TO_GALACTIC
Orient.DCM_MEME2000_TO_ECLIPJ2000
Orient.DCM_MEME2000_TO_B1950
Orient.DCM_MOON_PA421_TO_ME421
Orient.DCM_MOON_PA430_TO_ME421
Orient.DCM_MOON_PA430_TO_ME430
Orient.DCM_MOON_PA440_TO_ME421
```