# [Tempo](@id basic_tempo)

The `Tempo` module is thought to be a fast, efficient and precise time transformation library
capable to handle the different transformations needed in the astronomical and astrodynamical 
application of the `Basic` library.

## Overview

There are different ways to represent an epoch within `Basic`, depending on the specific 
application. This section is here to help you choose the proper time representation and to 
present the capabilities of the module in trasforming time between different representations.

First of all, there is a deep difference in the way time is tought in the everyday life and when dealing with space-related applications. Whenever we say _the 12:35 of the 1st of January 2023_, we are merging two concepts: the **calendar** (_1st January 2023_) and the **time representation** (_12:35_).  

### Calendars
A calendar is a system of organizing days. This is done by giving names to periods of time, typically days, weeks, months and years. A date is the designation of a single and specific day within such a system. The simplest calendar system just counts time periods from a reference date. This applies for the _Julian day_ or Unix Time.

##### Julian Calendar
The Julian day is the continuous count of days since the beginning of a Julian period. 
The Julian period is a chronological interval of 7980 years; year 1 of the Julian period was 
4713 BC (−4712). The Julian calendar year 2023 is year 6736 of the current Julian period. 
The next Julian Period begins in the year AD 3268.

##### Gregorian Calendar
The Gregorian calendar, introduced in 1582 by Pope Gregory XIII, was a modification of the Julian calendar. The main change made in the Gregorian calendar was the introduction of a new 
rule for *leap years*. In the Julian calendar, every fourth year was a leap year, regardless of whether or not it was a century year (i.e. a year ending in 00). However, the solar year is actually about 11 minutes shorter than 365.25 days, which means that over time, the Julian calendar began to drift away from the solar year. To correct this, the Gregorian calendar introduced a new rule that stated that years that are divisible by 100 are not leap years, unless they are also divisible by 400.

### Time & Scales

Calculations in any scientific discipline may involve precise time, but what 
sets astronomy apart is the number and variety of **time scales** that have to be used.
In fact, in astronomical applications the physical context of the “clock” matters,
whether it is on Earth, moving or stationary, or on a spacecraft.

```@raw html
<figure>
    <img src="https://gitlab.com/astronaut-tools/julia/Documentation/-/raw/390f98f53a0d35a3c0963dff8a5f608ff79304db/docs/src/assets/figures/enciclopedia/timescales.png" alt="Image" width="600" />
    <figcaption> Time Conversions - The difference in each timescale is shown with 
    respect to TAI. </figcaption>
</figure>
```

The most relevant time scales for these applications are:

- **UT1** (Universal Time 1): UT1 is a time scale based on the rotation of the Earth. 
    It is used to measure the positions of celestial objects relative to the Earth's 
    surface. UT1 is closely related to *Greenwich Mean Time (GMT)*, and the two time 
    scales are often used interchangeably.

- **TAI** (International Atomic Time): TAI is a time scale based on the average 
    frequency of a set of atomic clocks. It is used to measure the positions of 
    celestial objects relative to the Earth's surface.

- **TT** (Terrestrial Time): TT is a time scale based on the motion of celestial 
    objects around the solar system barycenter (the center of mass of the solar system). 
    It is used to measure the positions of celestial objects relative to the Earth's surface.

- **TDB** (Barycentric Dynamical Time): TDB is a time scale based on the motion of 
    celestial objects around the solar system barycenter (the center of mass of the 
    solar system). It is used to measure the positions of celestial objects relative 
    to the solar system barycenter.

- **TCB** (Barycentric Coordinate Time): TCB is a time scale based on the motion of 
    celestial objects around the solar system barycenter (the center of mass of the 
    solar system). It is used to measure the positions of celestial objects relative 
    to the solar system barycenter.

- **TCG** (Geocentric Coordinate Time): TCG is a time scale based on the rotation of 
    the Earth. It is used to measure the positions of celestial objects relative to the 
    Earth's surface.

- **Teph** (Ephemeris Time): Teph is a time scale based on the motion of celestial 
    objects around the solar system barycenter (the center of mass of the solar 
    system). It is used to measure the positions of celestial objects relative to 
    the solar system barycenter.

Of the seven time scales to be described here, one is atomic time (TAI), 
one is solar time (UT1), one is an atomic/solar hybrid (UTC) and four are 
dynamical times (TT, TCG, TCB, TDB). Other time scales of interest may also be the 
ones associated to the different positioning systems. In particular: **GPS** (Global 
Positioning System), **GLONASS** (Global Navigation Satellite System) and **GALILEO** 
(Global Navigation Satellite System) times could be defined as a constant offset with
respect to TAI.

## Time in `Basic`

Within `Basic`, the way in which time is represented in `Basic` is through the use of 
[`Epoch`](@ref)s. `Epoch`s are an efficient, differentiable and precise way to represent 
astronomical time. To parse an epoch object, two parameters shall be assigned:

- **Timescale**: This parameter determines the timescale that the epoch is based on. 
    For example, it can be set to UTC, TAI, TDB, or TCB. This allows the user to convert 
    the epoch between different timescales if necessary.

- **Origin**: This parameter determines the origin of the epoch, which is the point in time 
    from which the epoch is measured. This can be in the form of a Julian date, a
    Modified Julian date or any user-defined origin. 
    The origin can also be set to a specific event, such as J2000.0 or B1950.0.

By assigning these two parameters, [`Epoch`](@ref)s can be used to represent time in a precise 
manner, which is crucial for accurate timekeeping and coordination of events in a universe model.

## [API](@id basic_tempo_api)

### Types

```@meta
DocTestSetup = quote
    using Basic
end
```

```@autodocs
Modules = [Basic.Tempo]
Order = [:type]
```

### Functions

```@autodocs
Modules = [Basic.Tempo]
Order = [:function]
```

### Macros

```@autodocs
Modules = [Basic.Tempo]
Order = [:macro]
```

### Constants

```@autodocs
Modules = [Basic.Tempo]
Order = [:constant]
``` 