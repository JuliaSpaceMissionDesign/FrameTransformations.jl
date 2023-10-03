export add_axes_eclipj2000!, add_axes_meme2000!, add_axes_mememod!

"""
    add_axes_meme2000!(frames::FrameSystem{O,T}, axes::AbstractFrameAxes, parent::AbstractFrameAxes) where {O, T}

Add `axes` as a set of inertial axes representing the Mean Equator Mean Equinox of J2000 
to `frames`. 

!!! warning 
    The name (or the axes ID) of the parent set of axes must be `ICRF` (i.e., the International 
    Celestial Reference Frame, ID = 1), or the `ECLIPJ2000` (i.e., the Ecliptic Equinox of 7
    J2000, ID = 17), otherwise and error is thrown.

### Examples 
```julia-repl 
julia> @axes ICRF 1 InternationalCelestialReferenceFrame

julia> @axes MEME2000 22 MeanEquatorMeanEquinoxJ2000 

julia> FRAMES = FrameSystem{1, Float64}();

julia> add_axes_inertial!(FRAMES, ICRF)

julia> add_axes_meme2000!(FRAMES, MEME2000, ICRF)
```

### See also 
See also [`add_axes_inertial!`](@ref) and [`Orient.DCM_ICRF_TO_J2000_BIAS`](@ref)

----

    add_axes_meme2000!(frames::FrameSystem{O,T}, name::Symbol, parentid::Int) where {T,O}

Add a new inertial axes representing Mean Equator Mean Equinox of J2000 giving a `name` and 
a `parent`. The id is automatically assigned as [`Orient.AXESID_MEME2000`](@ref).

"""
function add_axes_meme2000!(
    frames::FrameSystem{O,T}, name::Symbol, parentid::Int
) where {T,O}

    if parentid == Orient.AXESID_ICRF
        dcm = Orient.DCM_ICRF_TO_J2000_BIAS
    elseif parentid == Orient.AXESID_ECLIPJ2000
        dcm = Orient.DCM_J2000_TO_ECLIPJ2000'
    else 
        throw(
            ArgumentError(
                "Mean Equator, Mean Equinox of J2000 (MEME2000) axes can only be defined " *
                "w.r.t. the International Celestial Reference Frame (ICRF).",
            ),
        )
    end
    return add_axes_inertial!(frames, name, Orient.AXESID_MEME2000; parentid=parentid, dcm=dcm)
end

function add_axes_meme2000!(
    frames::FrameSystem{O,T}, axes::AbstractFrameAxes, parent::AbstractFrameAxes
) where {T,O}
    pname = axes_name(parent)
    pid = axes_id(parent)

    if pname == :ICRF || pid == Orient.AXESID_ICRF
        dcm = Orient.DCM_ICRF_TO_J2000_BIAS

    elseif pname == :ECLIPJ2000 || pid == Orient.AXESID_ECLIPJ2000
        dcm = Orient.DCM_J2000_TO_ECLIPJ2000'
    else
        throw(
            ArgumentError(
                "Mean Equator, Mean Equinox of J2000 (MEME2000) axes can only be defined " *
                "w.r.t. the International Celestial Reference Frame (ICRF).",
            ),
        )
    end
    return add_axes_inertial!(frames, axes; parent=parent, dcm=dcm)
end

"""
    add_axes_eclipj2000!(frames, axes, parent::AbstractFrameAxes, iau_model::IAUModel=iau1980)
    
Add `axes` as a set of inertial axes representing the Ecliptic Equinox of J2000 (ECLIPJ2000)
to `frames`. The obliquity of the ecliptic is computed using the IAU Model `iau_model`.

The admissed `parent` set of axes are the following: 
- **ICRF**: for the International Celestial Reference Frame, with ID = 1
- **MEME2000**: the Mean Earth/Moon Ephemeris of J2000, with ID = 22

!!! warning 
    If the name (or the axes ID) of the parent set of `axes` is neither ICRF (ID = 1) nor 
    MEME2000 (ID = 22), an error is thrown. 

### Examples
```julia-repl 
julia> @axes ICRF 1 InternationalCelestialReferenceFrame

julia> @axes ECLIPJ2000 17 EclipticEquinoxJ2000 

julia> FRAMES = FrameSystem{1, Float64}();

julia> add_axes_inertial!(FRAMES, ICRF)

julia> add_axes_eclipj2000!(FRAMES, ECLIPJ2000, ICRF)
``` 

### See also 
See also [`add_axes_inertial!`](@ref) and [`Orient.DCM_ICRF_TO_J2000_BIAS`](@ref)

----

    add_axes_eclipj2000!(frames::FrameSystem{O,T}, name::Symbol, parentid::Int) where {T,O}

Add a new inertial axes representing the Ecliptic Equinox of J2000 giving a `name` and 
a `parent`. The axesid is automatically assigned as [`Orient.AXESID_ECLIPJ2000`](@ref).

"""
function add_axes_eclipj2000!(
    frames::FrameSystem{O,T}, name::Symbol, parentid::Int
) where {T,O}

    j2000_to_eclip = angle_to_dcm(orient_obliquity(iau_model, 0.0), :X)

    if parentid == Orient.AXESID_ICRF
        dcm = j2000_to_eclip * Orient.DCM_ICRF_TO_J2000_BIAS
    elseif parentid == Orient.AXESID_MEME2000
        dcm = j2000_to_eclip
    else 
        throw(
            ArgumentError(
                "Ecliptic Equinox of J2000 (ECLIPJ2000) axes could not be defined" *
                " w.r.t. $parentid axes. Only `ICRF` ($(Orient.AXESID_ICRF)) or" * 
                " `MEME2000`($(Orient.AXESID_MEME2000)) are accepted as parent axes.",
            ),
        )
    end
    return add_axes_inertial!(frames, name, Orient.AXESID_ECLIPJ2000; parentid=parentid, dcm=dcm)
end

function add_axes_eclipj2000!(
    frames::FrameSystem{O,T},
    axes::AbstractFrameAxes,
    parent::AbstractFrameAxes,
    iau_model::Orient.IAUModel=Orient.iau1980,
) where {T,O}
    pname = axes_name(parent)
    pid = axes_id(parent)

    # Compute ecliptic orientation using the specified obliquity model!
    j2000_to_eclip = angle_to_dcm(orient_obliquity(iau_model, 0.0), :X)

    if pname == :ICRF || pid == Orient.AXESID_ICRF
        dcm = j2000_to_eclip * Orient.DCM_ICRF_TO_J2000_BIAS
    elseif pname == :MEME2000 || pid == Orient.AXESID_MEME2000
        dcm = j2000_to_eclip
    else
        throw(
            ArgumentError(
                "Ecliptic Equinox of J2000 (ECLIPJ2000) axes could not be defined" *
                " w.r.t. $pname axes. Only `ICRF` or `MEME2000` are accepted as parent axes.",
            ),
        )
    end

    return add_axes_inertial!(frames, axes; parent=parent, dcm=dcm)
end

"""
    add_axes_mememod!(frames, axes, parent::AbstractFrameAxes, model::IAU2006Model=iau2006b)

Add `axes` as a set of projected axes representing the Mean of Date Ecliptic Equinox to 
`frames`. 

!!! note
    Despite this class of axes has a rotation matrix that depends on time, its derivatives 
    are assumed null.

!!! warning 
    The name of the `parent` set of axes must be the ICRF or have ID = 1 (i.e., the 
    International Celestial Reference Frame), otherwise an error is thrown. 

"""
function add_axes_mememod!(
    frames::FrameSystem, axes::AbstractFrameAxes, parent::AbstractFrameAxes
)
    if axes_name(parent) != :ICRF && axes_id(parent) != Orient.AXESID_ICRF
        throw(
            ArgumentError(
                "Mean Equator, Mean Equinox of date axes can only be defined " *
                "w.r.t. the International Celestial Reference Frame (ICRF)",
            ),
        )
    end

    return add_axes_projected!(frames, axes, parent, Orient.orient_rot3_icrf_to_mememod)
end
