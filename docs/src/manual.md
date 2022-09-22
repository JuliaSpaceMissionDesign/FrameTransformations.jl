# Overview

`Basic` is the most _basic_ tool, dedicated to space missions analysis and design, 
available in the `astronaut` context . 

## [Required features](@id basic_manual_req_features)

In the following table the desired list of features considered in the development of 
this package are reported, together with their development status:


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
    <th class="tg-fymr">FID</th>
    <th class="tg-fymr">Feature</th>
    <th class="tg-fymr">Module</th>
    <th class="tg-fymr">Status</th>
    <th class="tg-fymr">Stories</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-0pky">FIB1</td>
    <td class="tg-0pky">Read NASA JPL SPK type 1, 2, 21 files</td>
    <td class="tg-0pky">Ephemeris</td>
    <td class="tg-0pky">🟢</td>
    <td class="tg-0pky">
        <a href="https://gitlab.com/astronaut-tools/julia/core/Basic/-/issues/2">#2</a>
    </td>
  </tr>
  <tr>
    <td class="tg-0pky">FIB2</td>
    <td class="tg-0pky">Read NASA JPL constants files (tpc)</td>
    <td class="tg-0pky">Utils</td>
    <td class="tg-0pky">🟢</td>
    <td class="tg-0pky">
        <a href="https://gitlab.com/astronaut-tools/julia/core/Basic/-/issues/3">#3</a>
    </td>
  </tr>
  <tr>
    <td class="tg-0pky">FIB3</td>
    <td class="tg-0pky">Read NASA JPL frame kernels</td>
    <td class="tg-0pky"><i>TDB</i></td>
    <td class="tg-0pky">🔴</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">FIB4</td>
    <td class="tg-0pky">Write NASA JPL frame kernels</td>
    <td class="tg-0pky"><i>TDB</i></td>
    <td class="tg-0pky">🔴</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">FIB5</td>
    <td class="tg-0pky">
        Parse geometrical properties associated to celestial bodies. 
        Shall include at least:
        <ul>
            <li><b>dynamical</b>: gravitational parameter</li>
            <li><b>geometrical</b>: mean radius, polar radius, equatorial radius</li>
            <li><b>IAU constants</b>: right ascension, declination, polar motion 
                and their derivatives (at least > 1<sup>st</sup>)
            </li>
        </ul>
    </td>
    <td class="tg-0pky">Bodies, Orient</td>
    <td class="tg-0pky">🟡</td>
    <td class="tg-0pky">
        <a href="https://gitlab.com/astronaut-tools/julia/core/Basic/-/issues/3">#3</a>
    </td>
  </tr>
  <tr>
    <td class="tg-0pky">FIB6</td>
    <td class="tg-0pky">Parse epochs in different formats</td>
    <td class="tg-0pky">Tempo</td>
    <td class="tg-0pky">🟡</td>
    <td class="tg-0pky">
        <a href="https://gitlab.com/astronaut-tools/julia/core/Basic/-/issues/17">#17</a>
    </td>
  </tr>
  <tr>
    <td class="tg-0pky">FIB7</td>
    <td class="tg-0pky">
        Transform epochs between different timescales, supporting at least: TDB,
        TCB, TT, TAI, TCG, UTC, UT1, GPS
    </td>
    <td class="tg-0pky">Tempo</td>
    <td class="tg-0pky">🟡</td>
    <td class="tg-0pky">
        <a href="https://gitlab.com/astronaut-tools/julia/core/Basic/-/issues/15">#15</a>
    </td>
  </tr>
  <tr>
    <td class="tg-0pky">FIB8</td>
    <td class="tg-0pky">
        Construct IAU standard rotation matrices.
    </td>
    <td class="tg-0pky">Orient, Frames</td>
    <td class="tg-0pky">🟡</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">FIB9</td>
    <td class="tg-0pky">
        Read and parse IERS standard rotation matrices.
    </td>
    <td class="tg-0pky">Orient, Frames</td>
    <td class="tg-0pky">🔴</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">FIB10</td>
    <td class="tg-0pky">
        Support local orbital frames (i.e. synodic frame).
    </td>
    <td class="tg-0pky">Frames</td>
    <td class="tg-0pky">🔴</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">FIB11</td>
    <td class="tg-0pky">
        Support frozen frames.
    </td>
    <td class="tg-0pky">Frames</td>
    <td class="tg-0pky">🔴</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">FIB12</td>
    <td class="tg-0pky">
        Parse bodies and construct body graphs
    </td>
    <td class="tg-0pky">Bodies</td>
    <td class="tg-0pky">🟡</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">FIB13</td>
    <td class="tg-0pky">
        Perform frames transformations between any frame in the environment, with
        any center used as reference and any timescale used for epoch.
    </td>
    <td class="tg-0pky">Frames</td>
    <td class="tg-0pky">🔴</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">FIB14</td>
    <td class="tg-0pky">
        Read YAML, JSON environment configuration files.
    </td>
    <td class="tg-0pky">Universe, Utils</td>
    <td class="tg-0pky">🔴</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">FIB15</td>
    <td class="tg-0pky">
        Create universe autogenerated executable file from YAML, JSON.
    </td>
    <td class="tg-0pky">All</td>
    <td class="tg-0pky">🔴</td>
    <td class="tg-0pky"></td>
  </tr>
  <tr>
    <td class="tg-0pky">FIB16</td>
    <td class="tg-0pky">
        Read YAML, JSON, TEXT, JLD2 files with a simple interface to all.
    </td>
    <td class="tg-0pky">Utils</td>
    <td class="tg-0pky">🟡</td>
    <td class="tg-0pky"></td>
  </tr>
</tbody>
</table>
```

|    | Description                           |
|----|-------------------------------------- |
| 🔵 | Stable, deeply tested                 |
| 🟢 | Developed, working, partially tested  |
| 🟡 | In development                        |
| 🔴 | Development not started but planned   |
| ⚪ | Outdated/no more supported            |

# Design goals

The development of this package has been performed with the following design goals in mind:

1. **Efficiency**: being a fundamental package a particular attention has been 
    given to the execution time of the different routines as well as most/all of
    them have been optimised and deeply benchmarked.

2. **Extensibility**: attention has been given also to the definition of the 
    interfaces, which have been kept the most essential possible, in such a way 
    their extension can be performed very easily (also thanks to Julia language itself).

3. **Single Responsability Principle**: The different modules in this package 
    have been organized in such a way they are responsible of bringing only *one* 
    of the desired features. This results in the possibility to extend and maybe, 
    in future, detatch some modules to a different package.

# Universe configuration

The flow which is followed in this process is shown here:

```@raw html
<div class="mermaid">
flowchart LR

    classDef FixFont font-size:12px
    id1([config.yml]):::FixFont --> Read
    id2([config.json]):::FixFont --> Read 
    Read:::FixFont --> Schema:::FixFont
    Schema:::FixFont --> Read
    Read --> Parser:::FixFont --> id3((" ")) 

    subgraph pars [" "]
        id3((" ")) --> Utils:::FixFont
        id3((" ")) --> Constants:::FixFont 
        id3((" ")) --> Ephemeris:::FixFont 
        id3((" ")) --> Bodies:::FixFont 
        id3((" ")) --> Orient:::FixFont 
        id3((" ")) --> Frames:::FixFont 
    end

    pars --> id4([universe.jl]):::FixFont
</div>
```

Therefore first, the configuration input file is validated through the `Universe` schema,
then the `Parser` distribute to the different modules which parse a part of 
the an autogenerated file which shall be included then for the simulation 
environment to be loaded, i.e. via `include("universe.jl")`.