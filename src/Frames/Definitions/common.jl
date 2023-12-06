export add_axes_icrf!, add_axes_gcrf!, add_axes_meme2000!, add_axes_mod!, add_axes_tod!

"""
    add_axes_icrf!(frames::FrameSystem)

Add the International Celestial Reference Frame (ICRF) as the root axes of the frames graph.
The axes are automatically named `ICRF` and assigned the $(Orient.AXESID_ICRF) ID. 

### See also 
See also [`add_axes_inertial!`](@ref), [`add_axes_gcrf!`](@ref) and [`Orient.AXESID_ICRF`](@ref).
"""
@inline function add_axes_icrf!(frames::FrameSystem)

    if !isempty(frames_axes(frames))
        throw(ArgumentError("The ICRF can only be defined as a set of root axes."))
    end

    return add_axes_inertial!(frames, :ICRF, Orient.AXESID_ICRF)

end


"""
    add_axes_gcrf!(frames::FrameSystem)

Add the Geocentric Celestial Reference Frame (GCRF) to the frames graph. The axes are 
automatically named `GCRF` and assigned the $((Orient.AXESID_GCRF)) ID. These axes can only 
be defined as a set of root axes or as child of the ICRF (ID = $(Orient.AXESID_ICRF)).

### See also 
See also [`add_axes_inertial!`](@ref), [`add_axes_icrf!`](@ref) and [`Orient.AXESID_GCRF`](@ref).
"""
function add_axes_gcrf!(frames::FrameSystem)

    if has_axes(frames, Orient.AXESID_ICRF)
        # Add the GCRF as a child of the ICRF with an identity rotation 
        return add_axes_fixedoffset!(
            frames, :GCRF, Orient.AXESID_GCRF, Orient.AXESID_ICRF, DCM(1.0I)
        )

    elseif isempty(frames_axes(frames))
        # Add the GCRF as a root set of axes
        return add_axes_inertial!(frames, :GCRF, Orient.AXESID_GCRF)
        
    else 
        throw(
            ArgumentError(
                "The GCRF can only be defined with respect to the ICRF (ID =" * 
                " $(Orient.AXESID_ICRF)) or as a set of root axes."
            )
        )
    end
end


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
See also [`add_axes_inertial!`](@ref) and [`Orient.DCM_ICRF_TO_MEME2000`](@ref)
"""
@inline function add_axes_meme2000!(frames::FrameSystem, axes::AbstractFrameAxes, parent)
    return add_axes_meme2000!(frames, axes_name(axes), axes_alias(parent), axes_id(axes))
end

# Low-level function
function add_axes_meme2000!(
    frames::FrameSystem, name::Symbol, parentid::Int, axesid::Int=Orient.AXESID_MEME2000
)

    if parentid == Orient.AXESID_ICRF
        dcm = Orient.DCM_ICRF_TO_MEME2000
    elseif parentid == Orient.AXESID_ECLIPJ2000
        dcm = Orient.DCM_MEME2000_TO_ECLIPJ2000'
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
    add_axes_mod!(frames, axes::AbstractFrameAxes, parent)

Add `axes` as a set of projected axes representing the Mean Equator and Equinox of Date (MOD)
to `frames`. 

!!! note
    Despite this class of axes has a rotation matrix that depends on time, its derivatives 
    are assumed null.

!!! warning 
    The ID of the `parent` set of axes must be $(Orient.AXESID_ICRF) (ICRF), 
    otherwise an error is thrown. 

----

    add_axes_mod!(frames, name::Symbol, axesid::Int, parentid::Int)

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) type 
via the [`@axes`](@ref) macro.

### See also 
See also [`add_axes_projected!`](@ref) and [`Orient.orient_rot3_icrf_to_mod`](@ref)
"""
@inline function add_axes_mod!(frames::FrameSystem, axes::AbstractFrameAxes, parent)
    return add_axes_mod!(frames, axes_name(axes), axes_id(axes), axes_alias(parent))
end

# Low-level function
function add_axes_mod!(frames::FrameSystem, name::Symbol, axesid::Int, parentid::Int)

    if parentid != Orient.AXESID_ICRF
        throw(
            ArgumentError(
                "Mean of Date (MOD) axes can only be defined " *
                "w.r.t. the ICRF (ID = $(Orient.AXESID_ICRF)).",
            )
        )
    end

    return add_axes_projected!(
        frames, name, axesid, parentid, Orient.orient_rot3_icrf_to_mod
    )
end
