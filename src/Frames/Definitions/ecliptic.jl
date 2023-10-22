export add_axes_eclipj2000!, add_axes_meme2000!, add_axes_mememod!

"""
    add_axes_meme2000!(frames, axes::AbstractFrameAxes, parent)
    
Add `axes` as a set of inertial axes representing the Mean Equator Mean Equinox of J2000 
to `frames`. 

!!! warning 
    The the axes ID of the parent set of axes must be $(Orient.AXESID_ICRF) (ICRF) or 
    $(Orient.AXESID_ECLIPJ2000) (ECLIPJ2000) otherwise and error is thrown.

----

    add_axes_meme2000!(frames, name::Symbol, parentid::Int, axesid::Int = Orient.AXESID_MEME2000)

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) type 
via the [`@axes`](@ref) macro.

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
"""
@inline function add_axes_meme2000!(frames::FrameSystem, axes::AbstractFrameAxes, parent)
    return add_axes_meme2000!(frames, axes_name(axes), axes_alias(parent), axes_id(axes))
end

# Low-level function
function add_axes_meme2000!(
    frames::FrameSystem, name::Symbol, parentid::Int, axesid::Int=Orient.AXESID_MEME2000
)

    if parentid == Orient.AXESID_ICRF
        dcm = Orient.DCM_ICRF_TO_J2000_BIAS
    elseif parentid == Orient.AXESID_ECLIPJ2000
        dcm = Orient.DCM_J2000_TO_ECLIPJ2000'
    else 
        throw(
            ArgumentError(
                "Mean Equator, Mean Equinox of J2000 (MEME2000) axes can only be defined " *
                "w.r.t. the ICRF (ID = $(Orient.AXESID_ICRF)).",
            ),
        )
    end

    if axesid != Orient.AXESID_MEME2000
        @warn "$name is aliasing an ID that is not the standard MEME2000 ID" *
              " ($(Orient.AXESID_MEME2000))."
    end

    return add_axes_inertial!(frames, name, axesid; parentid = parentid, dcm = dcm)

end


"""
    add_axes_eclipj2000!(frames, axes::AbstractFrameAxes, parent, iau_model::IAUModel=iau1980)
    
Add `axes` as a set of inertial axes representing the Ecliptic Equinox of J2000 (ECLIPJ2000)
to `frames`. The obliquity of the ecliptic is computed using the IAU Model `iau_model`.

The admissed `parent` set of axes are the following: 
- **ICRF**: for the International Celestial Reference Frame, with ID = $(Orient.AXESID_ICRF)
- **MEME2000**: the Mean Earth/Moon Ephemeris of J2000, with ID = $(Orient.AXESID_MEME2000)

----

    add_axes_eclipj2000!(frames, name::Symbol, parentid::Int, iau_model::IAUModel=iau1980, 
        axesid::Int = Orient.AXESID_ECLIPJ2000)

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) type 
via the [`@axes`](@ref) macro.

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
"""
@inline function add_axes_eclipj2000!(
    frames::FrameSystem,
    axes::AbstractFrameAxes,
    parent,
    iau_model::Orient.IAUModel = Orient.iau1980,
)
    return add_axes_eclipj2000!(
        frames, axes_name(axes), axes_alias(parent), iau_model, axes_id(axes)
    )

end

# Low-level function
function add_axes_eclipj2000!(
    frames::FrameSystem, 
    name::Symbol, 
    parentid::Int, 
    iau_model::Orient.IAUModel = Orient.iau1980, 
    axesid::Int = Orient.AXESID_ECLIPJ2000
)

    # Compute the J2000 to ECLIPJ2000 rotationa ccording to the desired IAU model
    DCM_J2000_TO_ECLIPJ2000 = angle_to_dcm(orient_obliquity(iau_model, 0.0), :X)

    if parentid == Orient.AXESID_ICRF
        dcm = DCM_J2000_TO_ECLIPJ2000 * Orient.DCM_ICRF_TO_J2000_BIAS
    elseif parentid == Orient.AXESID_MEME2000
        dcm = DCM_J2000_TO_ECLIPJ2000
    else 
        throw(
            ArgumentError(
                "Ecliptic Equinox of J2000 (ECLIPJ2000) axes cannot be defined" *
                " w.r.t. $parentid axes. Only `ICRF` (ID = $(Orient.AXESID_ICRF)) or" * 
                " `MEME2000` (ID = $(Orient.AXESID_MEME2000)) are accepted as parent axes.",
            ),
        )
    end

    if axesid != Orient.AXESID_ECLIPJ2000
        @warn "$name is aliasing an ID that is not the standard ECLIPJ2000 ID" *
              " ($(Orient.AXESID_ECLIPJ2000))."
    end

    return add_axes_inertial!(
        frames, name, axesid; parentid = parentid, dcm = dcm)

end

"""
    add_axes_mememod!(frames, axes::AbstractFrameAxes, parent)

Add `axes` as a set of projected axes representing the Mean of Date Ecliptic Equinox to 
`frames`. 

!!! note
    Despite this class of axes has a rotation matrix that depends on time, its derivatives 
    are assumed null.

!!! warning 
    The ID of the `parent` set of axes must be $(Orient.AXESID_ICRF) (ICRF), 
    otherwise an error is thrown. 

----

    add_axes_eclipj2000!(frames, name::Symbol, axesid::Int, parentid::Int)

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) type 
via the [`@axes`](@ref) macro.

### See also 
See also [`add_axes_projected!`](@ref) and [`Orient.orient_rot3_icrf_to_mememod`](@ref)
"""
@inline function add_axes_mememod!(frames::FrameSystem, axes::AbstractFrameAxes, parent)
    return add_axes_mememod!(frames, axes_name(axes), axes_id(axes), axes_alias(parent))
end

# Low-level function
function add_axes_mememod!(frames::FrameSystem, name::Symbol, axesid::Int, parentid::Int)

    if parentid != Orient.AXESID_ICRF
        throw(
            ArgumentError(
                "Mean Equator, Mean Equinox of date axes can only be defined " *
                "w.r.t. the ICRF (ID = $(Orient.AXESID_ICRF)).",
            )
        )
    end

    return add_axes_projected!(
        frames, name, axesid, parentid, Orient.orient_rot3_icrf_to_mememod
    )
end
