# Welcome to Basic!

_Base tools for Mission Analysis, Navigation and Guidance._

Are you in search of fundamental routines to develop your own universe model? 
If so, `Basic` is the ideal starting point. The package is designed to provide users with 
the ability to create a customized, efficient, flexible, and extensible model for 
mission analysis and space mission design purposes. 

!!! compat 
    This package is still in development, some features may be modified 🙂

## Features 

- Convert between different time scales;
- Convert between different time representations;
- Read binary ephemeris files;
- Create a custom reference frame systems with user-defined points and axes.
- Transform states between different frames.  

All of this seamlessly integrated with `ForwardDiff.jl`.

## Mission

The development of this package has been performed with the following design goals in mind:

1. **Efficiency**: being a base package a particular attention has been 
    given to the execution time of the different routines as well as most/all of
    them have been optimised and deeply benchmarked.

2. **Extensibility**: attention has been given also to the definition of the 
    interfaces, which have been kept the most essential possible, in such a way 
    their extension can be performed very easily (also thanks to Julia language itself).

3. **Single Responsability Principle**: The different modules in this package 
    have been organized in such a way they are responsible of bringing only *one* 
    of the desired features. This results in the possibility to extend and maybe, 
    in future, detatch some modules to a different package.

4. **Automatic Differentiation**: seamless integration with `ForwardDiff.jl` is targetted 
    to fully exploit its power in higher-level packages constructed on top of `Basic`.


## [Requirements](@id basic_design_req)

The following table reports the desired list of features considered in the development of 
this package, together with their development status. 

```@raw html
<style type="text/css">
.tg  {border-collapse:collapse;border-spacing:0;}
.tg td{border-color:black;border-style:solid;border-width:1px;font-family:Arial, sans-serif;font-size:14px;
  overflow:hidden;padding:10px 5px;word-break:normal;}
.tg th{border-color:black;border-style:solid;border-width:1px;font-family:Arial, sans-serif;font-size:14px;
  font-weight:normal;overflow:hidden;padding:10px 5px;word-break:normal;}
.tg .tg-fymr{border-color:inherit;font-weight:bold;text-align:left;vertical-align:top}
.tg .tg-0pky{border-color:inherit;text-align:left;vertical-align:top}
</style>
<table class="tg">
<thead>
  <tr>
    <th class="tg-fymr">RID</th>
    <th class="tg-fymr">Feature</th>
    <th class="tg-fymr">Module</th>
    <th class="tg-fymr">Status</th>
    <th class="tg-fymr">Comment</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-0pky">REQ1</td>
    <td class="tg-0pky">Read binary JPL/INPOP ephemeris files. At least type 1, 2, 3, 21 shall be available.</td>
    <td class="tg-0pky">Ephemeris</td>
    <td class="tg-0pky">🔵</td>
    <td class="tg-0pky">
        The possibility to read binary SPK files is implemented through CALCEPH S/W. 
    </td>
  </tr>
  <tr>
    <td class="tg-0pky">REQ2</td>
    <td class="tg-0pky">Read JPL ASCII PCK files (constants).</td>
    <td class="tg-0pky">Utils</td>
    <td class="tg-0pky">🟢</td>
    <td class="tg-0pky">Fully implemented in Basic. Tested on NAIF's pck00010 and pck00011 files.</td>
  </tr>
  <tr>
    <td class="tg-0pky">REQ3</td>
    <td class="tg-0pky">Read NASA JPL frame kernels.</td>
    <td class="tg-0pky">Utils</td>
    <td class="tg-0pky">🔴</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">REQ4</td>
    <td class="tg-0pky">Write NASA JPL frame kernels</td>
    <td class="tg-0pky">Utils</td>
    <td class="tg-0pky">🔴</td>
    <td class="tg-0pky"></td>
  </tr>
    <tr>
    <td class="tg-0pky">REQ5</td>
    <td class="tg-0pky">
        Fetch automatically IERS EOP.
    </td>
    <td class="tg-0pky">Orient</td>
    <td class="tg-0pky">🟢</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">REQ6</td>
    <td class="tg-0pky">
        Fetch leapseconds files automatically.
    </td>
    <td class="tg-0pky">Orient</td>
    <td class="tg-0pky">🟢</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">REQ7</td>
    <td class="tg-0pky">Parse epochs in different formats. Shall include ISO and days or seconds since J2000.</td>
    <td class="tg-0pky">Tempo</td>
    <td class="tg-0pky">🔵</td>
    <td class="tg-0pky">
        Deeply tested against other S/W results.
    </td>
  </tr>
  <tr>
    <td class="tg-0pky">REQ8</td>
    <td class="tg-0pky">
        Transform epochs between different timescales, supporting at least: TDB,
        TCB, TT, TAI, TCG, UTC, UT1, GPS.
    </td>
    <td class="tg-0pky">Tempo</td>
    <td class="tg-0pky">🔵</td>
    <td class="tg-0pky">
        Deeply tested against ERFA. Extensible structure.
    </td>
  </tr>
  <tr>
    <td class="tg-0pky">REQ9</td>
    <td class="tg-0pky">
        Create IAU standard-based body-fixed rotation matrices.
    </td>
    <td class="tg-0pky">Orient</td>
    <td class="tg-0pky">🔵</td>
    <td class="tg-0pky">Deeply tested against NAIF's SPICE.</td>
  </tr>
  <tr>
    <td class="tg-0pky">REQ10</td>
    <td class="tg-0pky">
        Create ITRF (IERS-based) rotation matrices.
    </td>
    <td class="tg-0pky">Orient</td>
    <td class="tg-0pky">🔵</td>
    <td class="tg-0pky">Deeply tested against ERFA. Models with different precisions are available.</td>
  </tr>
  <tr>
    <td class="tg-0pky">REQ11</td>
    <td class="tg-0pky">
        Create a graph of custom points.
    </td>
    <td class="tg-0pky">Frames</td>
    <td class="tg-0pky">🔵</td>
    <td class="tg-0pky">Points could be defined in different ways (not only associated to ephemeris).</td>
  </tr>
  <tr>
    <td class="tg-0pky">REQ12</td>
    <td class="tg-0pky">
        Create a graph of custom axes.
    </td>
    <td class="tg-0pky">Frames</td>
    <td class="tg-0pky">🔵</td>
    <td class="tg-0pky">Axes models for MEME2000, ITRF, IAU models and many others are already implemented and tested within Basic.</td>
  </tr>
  <tr>
    <td class="tg-0pky">REQ13</td>
    <td class="tg-0pky">
        Get a state of a point relative to another in a custom frame.
    </td>
    <td class="tg-0pky">Frames</td>
    <td class="tg-0pky">🔵</td>
    <td class="tg-0pky">Deeply tested against NAIF's SPICE.</td>
  </tr>
  <tr>
    <td class="tg-0pky">REQ14</td>
    <td class="tg-0pky">
        Full compatibility with ForwardDiff.jl shall be assured.
    </td>
    <td class="tg-0pky">Basic</td>
    <td class="tg-0pky">🟡</td>
    <td class="tg-0pky">Partial compatibility is already available (see dedicated section for details)</td>
  </tr>
</tbody>
</table>
```

Where the color legend used is as follows:

|    | Description                           |
|----|-------------------------------------- |
| 🔵 | Stable, deeply tested                 |
| 🟢 | Developed, working, partially tested  |
| 🟡 | In development                        |
| 🔴 | Development not started but planned   |
| ⚪ | Outdated/no more supported            |

