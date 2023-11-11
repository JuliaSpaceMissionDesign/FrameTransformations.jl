# [Orient's Low-level API](@id low_orient_api)

Although this routines are not meant to be used outside of the package, they are here documented to aid future developments or to help users that require specific computations.

## Outdated IAU Models

```@docs
Orient.iau1980
```

## Fundamental Arguments 

```@docs 
Orient.FundamentalArguments
Orient.LuniSolarArguments
Orient.PlanetaryArguments
Orient.fa_mano_moon
Orient.fa_mano_sun
Orient.fa_mlon_moon
Orient.fa_mlat_moon
Orient.fa_melo_moon
Orient.fa_precession
Orient.fa_mlon_mercury
Orient.fa_mlon_venus
Orient.fa_mlon_earth
Orient.fa_mlon_mars
Orient.fa_mlon_jupiter
Orient.fa_mlon_saturn
Orient.fa_mlon_uranus
Orient.fa_mlon_neptune
```

## Precession

```@docs
Orient.frame_bias
Orient.fw_angles
Orient.fw_matrix
Orient.precession_angles
Orient.precession_rate
```

## IERS Routines

```@docs 
Orient.bpn2xy
Orient.cip_coords
Orient.cio_locator
Orient.cip_motion
Orient.earth_rotation_angle
Orient.earth_rotation_rate 
Orient.era_rotm
Orient.fw2xy
Orient.tio_locator
Orient.polar_motion
Orient.xys2m
```

## EOP Data

```@docs 
Orient.EOPData
Orient.EOPInterpolator
Orient.prepare_eop 
Orient.read_eop 
Orient.init_eop 
Orient.IERS_EOP
Orient.offset_utc2ut1
```

## Transformations

```@docs 
Orient.orient_rot3_icrf_to_mememod
```