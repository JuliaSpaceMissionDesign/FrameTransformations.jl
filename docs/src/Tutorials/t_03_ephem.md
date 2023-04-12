# [Ephemeris loading and manipulations](@id tutorial_03_ephem)

Differently from the other modules, `Ephemeris` is practically an interface-only package.
In this tutorial the most common use-case as well as the extensions capabilites are presented.


```julia
using Basic
```

## JPL/INPOP ephemeris handling

At the time being, the capabilities of [Calceph](https://www.imcce.fr/recherche/equipes/asd/calceph/)
C library is exploited via the [CALCEPH.jl](https://github.com/JuliaAstro/CALCEPH.jl) Julia
interface are used to access the binary planetary ephemeris files, such INPOPxx, JPL DExxx 
and SPICE ephemeris files.

!!! warning
    This introduces an inherent limitation, due to the fact that it is **not** possible to perform
    autodiff across C code called by Julia. Therefore, the general ephemeris handling may change
    in future, whenever a pure Julia ephemeris reader will be available.

Therefore, let's have a look to some use-cases in which this interface may be used.

First of all, the necessary kernels to follow-on the example shall be downloaded:
- [DE440](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp)
- [MOON\_PA\_DE440](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de440_200625.bpc)

Then, they can be loaded to the workspace exploiting the `CalcephProvider` interface:


```julia
eph = load(
    CalcephProvider, 
    [
        "/home/andrea/Documents/Kernels/pck/moon_pa_de440_200625.bpc", 
        "/home/andrea/Documents/Kernels/spk/de440.bsp"
    ]
)
```


    CalcephProvider(CALCEPH.Ephem(Ptr{Nothing} @0x0000000007049610))



```julia
# Get position records
ephem_position_records(eph)
```


    14-element Vector{CALCEPH.PositionRecord}:
     CALCEPH.PositionRecord(1, 0, 2.2871845e6, 2.6889765e6, 1)
     CALCEPH.PositionRecord(2, 0, 2.2871845e6, 2.6889765e6, 1)
     CALCEPH.PositionRecord(3, 0, 2.2871845e6, 2.6889765e6, 1)
     CALCEPH.PositionRecord(4, 0, 2.2871845e6, 2.6889765e6, 1)
     CALCEPH.PositionRecord(5, 0, 2.2871845e6, 2.6889765e6, 1)
     CALCEPH.PositionRecord(6, 0, 2.2871845e6, 2.6889765e6, 1)
     CALCEPH.PositionRecord(7, 0, 2.2871845e6, 2.6889765e6, 1)
     CALCEPH.PositionRecord(8, 0, 2.2871845e6, 2.6889765e6, 1)
     CALCEPH.PositionRecord(9, 0, 2.2871845e6, 2.6889765e6, 1)
     CALCEPH.PositionRecord(10, 0, 2.2871845e6, 2.6889765e6, 1)
     CALCEPH.PositionRecord(301, 3, 2.2871845e6, 2.6889765e6, 1)
     CALCEPH.PositionRecord(399, 3, 2.2871845e6, 2.6889765e6, 1)
     CALCEPH.PositionRecord(199, 1, 2.2871845e6, 2.6889765e6, 1)
     CALCEPH.PositionRecord(299, 2, 2.2871845e6, 2.6889765e6, 1)



```julia
# Get orientation records
ephem_orient_records(eph)
```


    2-element Vector{CALCEPH.OrientationRecord}:
     CALCEPH.OrientationRecord(31008, 2.2871845e6, 2.6071845e6, 1)
     CALCEPH.OrientationRecord(31008, 2.6071845e6, 2.6889765e6, 1)


Therefore, by means of the `ephem_load` method, we are capable to load multiple ephemeris files
with different extensions all in a single step. Note that the ephemeris loaded may be both
positions or orientations as they are both handled by `CalcephProvider`.

Now, once loaded, we can use the `eph` object to compute quantities as required. For example,
to compute the position of the Moon w.r.t. Saturn we may simply call `ephem_compute_order!`,
which is an in-place method for _arbitrary_ order ephemeris evaluation:


```julia
v = zeros(6)
ephem_compute!(v, eph, Basic.Tempo.DJ2000, 0.0, 301, 6, 0)
v[1:3]
```


    3-element Vector{Float64}:
     -9.85175766281782e8
     -7.912249575768421e8
     -2.828202495894639e8


If we want to compute the velocities then we just need to change the order to `1`, etc...


```julia
v = zeros(6)
ephem_compute!(v, eph, Basic.Tempo.DJ2000, 0.0, 301, 6, 1)
v
```


    6-element Vector{Float64}:
      -9.85175766281782e8
      -7.912249575768421e8
      -2.828202495894639e8
     -21.7187066847619
     -11.793316314234914
      -5.319653062269718


Note that, to use this function, the time shall be written as a two-part julian date, for
precision purposes.


The same approach may be used for orientations, where `ephem_orient_order!` can be used in a
similar way:


```julia
v = zeros(6)
ephem_orient!(v, eph, Basic.Tempo.DJ2000, 0.0, 31008, 0)
v[1:3]
```


    3-element Vector{Float64}:
       -0.054147058087019556
        0.42485546099906896
     2564.2582727215126


## Creating a custom ephemeris provider

As said, the `Ephemeris` package is more an interface than an actual implementation of an
ephemeris provider. Therefore, there is the possibility for the user to extend the 
functionalities of `Basic` with a custom ephemeris provider.

To that purpose, the following steps shall be performed:

1. Define a subtype of `AbstractEphemerisProvider`, creating the concrete implementation.
2. Define an overload of `ephem_load` to load data for the new ephemeris provider.
3. Define the required overloads of the abstract methods present in `Ephemeris\abstract.jl` (you don't have to define all of them, just the ones required for your specific application).

If the application includes the use of the new ephemeris provider within a `FrameSystem` object,
then additional overloads should be created to ease the `FrameSystem` parsing.
