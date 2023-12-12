export add_axes_eclipj2000!

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
See also [`add_axes_inertial!`](@ref) and [`Orient.DCM_ICRF_TO_MEME2000`](@ref)
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
    DCM_MEME2000_TO_ECLIPJ2000 = angle_to_dcm(orient_obliquity(iau_model, 0.0), :X)

    if parentid == Orient.AXESID_ICRF
        dcm = DCM_MEME2000_TO_ECLIPJ2000 * Orient.DCM_ICRF_TO_MEME2000
    elseif parentid == Orient.AXESID_MEME2000
        dcm = DCM_MEME2000_TO_ECLIPJ2000
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
