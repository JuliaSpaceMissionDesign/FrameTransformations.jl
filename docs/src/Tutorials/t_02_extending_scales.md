# [Timescales graphs and extensions](@id tutorial_02_scales)

In `Basic`, timescales are connected each other via a **directed graph**. Thanks to the 
structure of the `Tempo` module, it is possible to either extend the current graph of 
scales or create a completely custom one. Both of this possibilities are the subject of 
this tutorial.


```julia
using Basic.Tempo
```

## Create a timescales graph

To create a computational directed graph to handle timescales, `Basic` provides the `TimeSystem`
type. Therefore, let us define a new time transformation system called `TIMETRANSF`:


```julia
const TIMETRANSF = TimeSystem{Float64}()
```


    TimeSystem{Float64}(Basic.MappedGraphs.MappedNodeGraph{Basic.Tempo.TimeScaleNode{Float64}, Graphs.SimpleGraphs.SimpleDiGraph{Int64}}(Graphs.SimpleGraphs.SimpleDiGraph{Int64}(0, Vector{Int64}[], Vector{Int64}[]), Dict{Int64, Int64}(), Basic.Tempo.TimeScaleNode{Float64}[], Dict{Int64, Dict{Int64, Vector{Int64}}}(), Dict{Int64, Dict{Int64, Int64}}()))


This object contains a graph and the properties associated to the new time-system defined in
`TIMETRANSF`. Note that the computational graph at the moment is empty, thus, we need to 
manually populate it with the new transformations.

## Create a new timescale

In order to insert a new timescale to the graph, a new timescale type alias shall be defined.
This can be easily done via the macro `@timescale`. This step requires 3 elements:
- The timescale acronym (user-defined).
- The timescale index (it is an `Int` used to uniquely represent the timescale).
- The timescale fullname.


```julia
@timescale DTS 1 DefaultTimeScale
```

## Register the new timescale

Once, created, the new timescale is ready to be registered. If it is the first scale registered
in the computational graph, than, nothing else than the type alias is needed and the 
registration can be performed as follows:


```julia
add_timescale(TIMETRANSF, DTS)
```


```julia
TIMETRANSF.scales.nodes
```


    1-element Vector{Basic.Tempo.TimeScaleNode{Float64}}:
     TimeScaleNode{Float64}(name=DTS, id=1)



Instead, in case the timescale is linked to a parent one, an offset function shall be defined.
Remember that the computational graph is **directed**, i.e. the transformation to go back and 
forth to the parent shall be defined if two-way transformations are desired.

In this example, assume we want to register timescale `NTSA` and a timescale `NTSB`. 
`NTSA` has `DTS` as parent and a constant offset of 1 second. `NTSB` has `NTSA` has parent 
and a linear offset with slope of 1/86400.

Then, first create the new scales:


```julia
@timescale NTSA 2 NewTimeScaleA
@timescale NTSB 3 NewTimeScaleB
```

Now, let us define the offset functions for `NTSA`:


```julia
const OFFSET_DTS_TO_NTSA = 1.0
@inline offset_dts2ntsa(sec::Number) = OFFSET_DTS_TO_NTSA
@inline offset_ntsa2dts(sec::Number) = -OFFSET_DTS_TO_NTSA;
```

We can now register `NTSA` to the computational graph using the `add_timescale` method:


```julia
add_timescale(TIMETRANSF, NTSA, offset_dts2ntsa, parent=DTS, ftp=offset_ntsa2dts)
```

Now, if we have a look to the computational graph, we'll se that `NTSA` is registered:


```julia
TIMETRANSF.scales.nodes
```


    2-element Vector{Basic.Tempo.TimeScaleNode{Float64}}:
     TimeScaleNode{Float64}(name=DTS, id=1)
    
     TimeScaleNode{Float64}(name=NTSA, id=2, parent=1)



As well as, since we have registered both direct and inverse transformations, there is the 
possibility to transform back and forth from `NTSA` to `DTS`. We can easily see this looking
at the `paths` contained in the computational graph. Here the timescale are represented by means
of the type-alias unique integer assigned during the creation of the new type. 


```julia
TIMETRANSF.scales.paths
```


    Dict{Int64, Dict{Int64, Vector{Int64}}} with 2 entries:
      2 => Dict(1=>[2, 1])
      1 => Dict(2=>[1, 2])


If now we create a `DTS` epoch, it is possible to use the custom time transformation system
to convert to `NTSA`:


```julia
# Create the new epoch 
# IMPORTANT: only J2000 seconds Epoch parser works with custom timescales.
e = Epoch(0.0, DTS)
```


    2000-01-01T12:00:00.000 DTS



```julia
# Call `convert` using the custom time transformation system 
convert(NTSA, e, system=TIMETRANSF)
```


    2000-01-01T12:00:01.000 NTSA


!!! note
    The `system` is an optimal output if the `Basic` time transformation system is used.

To conclude the example, `NTSB` is has to be inserted. Let's assume that only the transformation
`NTSA -> NTSB` can be constructed. Then: 


```julia
# Create the linear offset function
offset_ntsa2ntsb(sec::Number) = sec/86400.0
# Register the timescale to the computational graph
add_timescale(TIMETRANSF, NTSB, offset_ntsa2ntsb, parent=NTSA)
```

Now, let's have a look to the nodes in the graph:


```julia
TIMETRANSF.scales.nodes
```


    3-element Vector{Basic.Tempo.TimeScaleNode{Float64}}:
     TimeScaleNode{Float64}(name=DTS, id=1)
    
     TimeScaleNode{Float64}(name=NTSA, id=2, parent=1)
    
     TimeScaleNode{Float64}(name=NTSB, id=3, parent=2)



Where the new timescale has been registered with the alias `3`. Note however, that from `3` no
transformations are available:


```julia
TIMETRANSF.scales.paths
```


    Dict{Int64, Dict{Int64, Vector{Int64}}} with 3 entries:
      2 => Dict(3=>[2, 3], 1=>[2, 1])
      3 => Dict(2=>[], 1=>[])
      1 => Dict(2=>[1, 2], 3=>[1, 2, 3])


To conclude, let's test the new time transformation system. Let's take the previous `Epoch` 
translate forward of 2 days and transform to `NTSA` and `NTSB`. We should obtain a translation
of `1 sec` and `3 sec` respectively:


```julia
# Translate the epoch
e += 2*86400
```


    2000-01-03T12:00:00.000 DTS



```julia
# Convert to `NTSA`
ea = convert(NTSA, e, system=TIMETRANSF)
```


    2000-01-03T12:00:01.000 NTSA



```julia
# Convert to `NTSB`
eb = convert(NTSB, e, system=TIMETRANSF)
```


    2000-01-03T12:00:03.000 NTSB

