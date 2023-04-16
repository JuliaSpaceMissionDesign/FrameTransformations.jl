# Contributing

Although `Basic` already implements may of its required functionalities, there is still 
a lot to do to finish it. 


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
    <th class="tg-fymr">ID</th>
    <th class="tg-fymr">Feature to develop</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td class="tg-0pky">AD1</td>
    <td class="tg-0pky">
        <b>Ephemeris</b>: the actual default ephemeris provider for Basic is Calceph. 
        While it is a very versatile library, there is a inherent issue in that: the 
        impossibility to perform AD over it. Therefore, within the context of development 
        of Basic a (full) julia reader for JPL/INPOP ephemeris is foreseen. This package
        is actually at the early stage of the development at the moment, thus any 
        contribution would be highly appreciated.

        This is being performed in a separate package called `Ephemeris` within the JSMD
        environment.
    </td>
  </tr>
    <tr>
    <td class="tg-0pky">AD2</td>
    <td class="tg-0pky">
        <b>Caching</b>: if you have a deep look into the FrameSystem structure, you'll 
        see that some frames have a dedicated cache. This goes in contrast with the use 
        of ForwardDiff, since it would require a Dual's compatible cache: there is, then,
        the need of replacing those caches with a PreallocationTools.DiffCache or something
        similar.
    </td>
  </tr>
  </tr>
    <tr>
    <td class="tg-0pky">AD3</td>
    <td class="tg-0pky">
        <b>FunctionWrappers</b>: if you have a deep look into the FrameSystem structure, you'll 
        see that the actual transformations are stored by means of FunctionWrapper. Again,
        this goes against the AD integration but was choosen for efficiency reasons.
        In this case there is the need to introduct FunctionWrappersWrappers dependency to
        solve the issue. 
    </td>
  </tr>
</tbody>
</table>
```



Write down an [email](mailto:andrea.pasquale@outlook.it) if your are interested in 
contributing to one of those, we'll be for welcoming you with some additional guidelines. 

## Guidelines

This section details the some of the guidelines that should be followed when contributing 
to this package.

Since the package is not consolidated yet, the most straight forward way to contribute is 
by creating a branch of the `master` and develop there the desired feature. To this step, 
follows the creation of a PR that will be accepted only if the new feature are **validated** 
by tests.