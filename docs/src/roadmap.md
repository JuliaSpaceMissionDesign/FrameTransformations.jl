
# Development Roadmap

The following table reports the desired list of features considered in the development of 
this package, together with their development status. The following color legend is used: 

|    | Description                           |
|----|-------------------------------------- |
| 🔵 | Stable, deeply tested                 |
| 🟢 | Developed, working, partially tested  |
| 🟡 | In development                        |
| 🔴 | Development not started but planned   |
| ⚪ | Outdated/no more supported            |


## Frames 

All the initially desired features of Frames have been implemented. We are currently working 
to avoid undesired allocations when using dual numbers.

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
    <th class="tg-fymr" style="width:40%">Feature</th>
    <th class="tg-fymr" style="width:10%">Status</th>
    <th class="tg-fymr" style="width:50%">Comment</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-0pky">
        Create a graph of custom points.
    </td>
    <td class="tg-0pky">🔵</td>
    <td class="tg-0pky">Points could be defined in different ways (not only associated to ephemeris).</td>
  </tr>
  <tr>
    <td class="tg-0pky">
        Create a graph of custom axes.
    </td>
    <td class="tg-0pky">🔵</td>
    <td class="tg-0pky">Axes models for MEME2000, ITRF, IAU models and many others are already implemented and tested within FrameTransformationsFrameTransformations.</td>
  </tr>
  <tr>
    <td class="tg-0pky">
        Get a state of a point relative to another in a custom frame.
    </td>
    <td class="tg-0pky">🔵</td>
    <td class="tg-0pky">Deeply tested against NAIF's SPICE.</td>
  </tr>
    <tr>
    <td class="tg-0pky">
        Full compatibility with ForwardDiff.jl shall be assured.
    </td>
    <td class="tg-0pky">🟡</td>
    <td class="tg-0pky">Partial compatibility is already available. We are resolving the 
    last issues related to undesired allocations when using dual numbers. </td>
  </tr>
</tbody>
</table>
```


## Orient

All the features of the Orient submodule have been developed and tested. No major updates
are currently foreseen.

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
    <th class="tg-fymr" style="width:40%">Feature</th>
    <th class="tg-fymr" style="width:10%">Status</th>
    <th class="tg-fymr" style="width:50%">Comment</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-0pky">
        Create IAU standard-based body-fixed rotation matrices.
    </td>
    <td class="tg-0pky">🔵</td>
    <td class="tg-0pky">Deeply tested against NAIF's SPICE.</td>
  </tr>
  <tr>
    <td class="tg-0pky">
        Create ITRF (IERS-based) rotation matrices.
    </td>
    <td class="tg-0pky">🔵</td>
    <td class="tg-0pky">Deeply tested against ERFA. Models with different precisions are available.</td>
  </tr>
  <tr>
    <td class="tg-0pky">
        Fetch leapseconds files automatically.
    </td>
    <td class="tg-0pky">🟢</td>
    <td class="tg-0pky"></td>
  </tr>
</tbody>
</table>
```
