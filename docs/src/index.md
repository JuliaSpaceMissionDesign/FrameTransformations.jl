# Welcome to Basic!

Are you in search of fundamental routines to develop your own universe model? If so, `Basic` is the ideal starting point. The package is designed to provide users with the ability to create a customized, efficient, flexible, and extensible multibody model for mission analysis and space mission design purposes. 

Specifically, `Basic` offers a wide range of functionality that can be used to create an efficient universe model. One of the key features of the package is its ability to transform time between different time scales and representations. This includes the ability to **convert between various time systems** such as UTC, TAI, TDB, and TCB. Additionally, `Basic` can also convert between different time representations such as Julian dates, Modified Julian dates, and calendar dates. This functionality is crucial for accurate timekeeping and coordination of events in a universe model.

Another important feature of Basic is its ability to **read JPL/INPOP ephemeris files**. These files contain highly accurate data on the positions and velocities of celestial objects such as planets and moons. By providing the ability to read these files, `Basic` allows users to incorporate this data into their universe model and ensure a high degree of accuracy.

`Basic` also provides functionality for **transforming states between different frames**. This includes the ability to convert between inertial and non-inertial frames of reference, as well as the ability to transform position, velocity, acceleration and jerk between different coordinate systems. This is important for modeling complex dynamics and interactions in a universe model.

!!! compat 
    This package is still in development, some features may be not available 
    or modified in the future 🙂

## Index

### Manual

```@contents
Pages = ["manual.md"]
Depth = 2
```

### Modules

```@contents
Pages = ["Modules/time.md", 
         "Modules/orient.md",
         "Modules/ephem.md",
         "Modules/utils.md"]
Depth = 2
```