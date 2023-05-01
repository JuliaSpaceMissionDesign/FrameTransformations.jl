# FrameSystem Benchmarks

```@setup benchmark_case
using Basic
using BenchmarkTools

# Create the new axes 
@axes ICRF 1 InternationalCelestialReferenceFrame

# Create frame system 
eph = load(
    CalcephProvider, 
    [
        "/home/andrea/Documents/Kernels/spk/de440.bsp", 
        "/home/andrea/Documents/Kernels/pck/moon_pa_de440_200625.bpc"
    ]
);
iau = Basic.load(TPC("/home/andrea/Documents/Kernels/pck/pck00010.tpc"));
fs = FrameSystem{3, Float64}(eph);

# Register the new axes in the graph as root axes
add_axes_inertial!(fs, ICRF)

# Add all points to FrameSystem
@point SSB 0
@point MEB 1 
@point VEB 2 
@point EMB 3 
@point MAB 4 
@point JUB 5
@point SAB 6
@point URB 7 
@point NEB 8 
@point SUN 10
@point MON 301
@point EAR 399
@point MER 199
@point VEN 299

add_point_root!(fs, SSB, ICRF)
add_point_ephemeris!(fs, MEB)
add_point_ephemeris!(fs, VEB)
add_point_ephemeris!(fs, EMB)
add_point_ephemeris!(fs, MAB)
add_point_ephemeris!(fs, JUB)
add_point_ephemeris!(fs, SAB)
add_point_ephemeris!(fs, URB)
add_point_ephemeris!(fs, NEB)
add_point_ephemeris!(fs, SUN)
add_point_ephemeris!(fs, MON)
add_point_ephemeris!(fs, EAR)
add_point_ephemeris!(fs, MER)
add_point_ephemeris!(fs, VEN)

# Add axes
@axes IAU_EARTH 3 
@axes IAU_MOON 4
@axes LME2000 5
@axes ITRF 6 
@axes MOONPA_DE440 31008 

add_axes_bcrtod!(fs, iau, EAR, IAU_EARTH, ICRF);
add_axes_bcrtod!(fs, iau, MON, IAU_MOON, ICRF); 
add_axes_bci2000!(fs, iau, MON, LME2000, ICRF);
add_axes_itrf!(fs, ITRF, ICRF) # default IAUModel is iau2006b
add_axes_pa440!(fs, MOONPA_DE440, ICRF)
```

## Case I: `CalcephProvider` read

**Setup**: Load _de440.bsp_ kernel and construct the frame system with all the points within 
the file. Extract the points and their centers from the position records.

- **1.1**: Direct (two-nodes) ephemeris file read, ephemeris axes:

```@repl benchmark_case
@benchmark vector9($fs, $SSB, $EMB, $ICRF, $(rand(0.0:1e8)))
```

- **1.2**: Three-nodes ephemeris file read, ephemeris axes:

```@repl benchmark_case
@benchmark vector9($fs, $SSB, $EAR, $ICRF, $(rand(0.0:1e8)))
```

- **1.3**: Four-nodes ephemeris file read, ephemeris axes:

```@repl benchmark_case
@benchmark vector9($fs, $EAR, $JUB, $ICRF, $(rand(0.0:1e8)))
```

## Case II: `CalcephProvider` read and transform

- **2.1**: Direct (two-nodes) ephemeris file read, single axes transformation:

```@repl benchmark_case
@benchmark vector9($fs, $SSB, $EMB, $IAU_EARTH, $(rand(0.0:1e8)))
```

- **2.2**: Three-nodes ephemeris file read, single axes transformation:

```@repl benchmark_case
@benchmark vector9($fs, $SSB, $EAR, $IAU_EARTH, $(rand(0.0:1e8)))
```

- **2.3**: Four-nodes ephemeris file read, single axes transformation:

```@repl benchmark_case
@benchmark vector9($fs, $JUB, $EAR, $IAU_EARTH, $(rand(0.0:1e8)))
```

- **2.4**: Four-nodes ephemeris file read, single ephemeris axes transformation:

```@repl benchmark_case
@benchmark vector9($fs, $JUB, $EAR, $MOONPA_DE440, $(rand(0.0:1e8)))
```

## Case III: Rotation chains

- **3.1**: Single axes transformation (two-nodes), time-dependant axes:

```@repl benchmark_case
@benchmark rotation9($fs, $ICRF, $IAU_EARTH, $(rand(0.0:1e8)))
```

- **3.2**: Three-nodes axes transformation, time-dependant axes:

```@repl benchmark_case
@benchmark rotation9($fs, $IAU_MOON, $IAU_EARTH, $(rand(0.0:1e8)))
```

- **3.3**: Three-nodes axes transformation, high-precision time-dependant axes:

```@repl benchmark_case
@benchmark rotation9($fs, $IAU_MOON, $ITRF, $(rand(0.0:1e8)))
```

- **3.4**: Three-nodes axes transformation, high-precision time-dependant & ephemeris axes:

```@repl benchmark_case
@benchmark rotation9($fs, $MOONPA_DE440, $ITRF, $(rand(0.0:1e8)))
```