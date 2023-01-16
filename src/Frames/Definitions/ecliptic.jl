export add_axes_eclipj2000!, add_axes_meme2000!, add_axes_mememod!


"""
    add_axes_meme2000!(frames::, axes, parent::AbstractFrameAxes)

Add `axes` as a set of inertial axes representing the Mean Equator Mean Equinox of J2000 
to `frames`. 

!!! warning 
    The `parent` set of axes must be the International Celestial Reference Frame (ICRF), 
    otherwise and error is thrown.

### Examples 
```jldoctest 
julia> @axes ICRF 1 InternationalCelestialReferenceFrame

julia> @axes MEME2000 2 MeanEquatorMeanEquinoxJ2000 

julia> add_axes_inertial!(FRAMES, ICRF)

julia> add_axes_meme2000!(FRAMES, MEME2000, ICRF)
```

### See also 
See also [`add_axes_inertial!`](@ref) and [`Orient.DCM_ICRF2J2000_BIAS`](@ref)
"""
function add_axes_meme2000!(frames::FrameSystem{O, T}, axes::AbstractFrameAxes, 
            parent::AbstractFrameAxes) where {T, O}

    if axes_name(parent) != :ICRF 
        throw(
            ErrorException(
                "Mean Equator, Mean Equinox of J2000 (MEME2000) axes can only be defined "*
                "w.r.t the International Celestial Reference Frame (ICRF)."
            )
        )
    end

    add_axes_inertial!(frames, axes; parent=parent, dcm=Orient.DCM_ICRF2J2000_BIAS)
end


"""
    add_axes_eclipj2000!(frames, axes, parent::AbstractFrameAxes)
    
Add `axes` as a set of inertial axes representing the Ecliptic Equinox of J2000 (ECLIPJ2000)
to `frames`.

The `parent` set of axes can be either the International Celestial Reference Frame (ICRF) 
or the Mean Earth/Moon Ephemeris of 2000 (MEME2000). If the `parent` set of axes is ICRF, 
the orientation of the ECLIPJ2000 axes is defined using a rotation matrix called 
[`DCM_J20002ECLIPJ2000`](@ref) and a bias matrix called [`DCM_ICRF2J2000_BIAS`](@ref). 
If the parent set of axes is MEME2000, the orientation of the ECLIPJ2000 axes is defined using 
only the [`DCM_J20002ECLIPJ2000`](@ref) rotation matrix.

!!! warning 
    If the parent set of `axes` is neither ICRF nor MEME2000, an error is thrown. 

### Examples
```jldoctest 
julia> @axes ICRF 1 InternationalCelestialReferenceFrame

julia> @axes ECLIPJ2000 17 EclipticEquinoxJ2000 

julia> add_axes_inertial!(FRAMES, ICRF)

julia> add_axes_eclipj2000!(FRAMES, MEME2000, ICRF)
``` 

### See also 
See also [`add_axes_inertial!`](@ref), [`Orient.DCM_J20002ECLIPJ2000`](@ref) and 
[`Orient.DCM_ICRF2J2000_BIAS`](@ref)
"""
function add_axes_eclipj2000!(frames::FrameSystem{O, T}, axes::AbstractFrameAxes,
            parent::AbstractFrameAxes) where {T, O}

    pname = axes_name(parent)

    if pname == :ICRF 
        dcm = Orient.DCM_J20002ECLIPJ2000 * Orient.DCM_ICRF2J2000_BIAS 
    elseif pname == :MEME2000 
        dcm = Orient.DCM_J20002ECLIPJ2000
    else
        throw(
            ErrorException("Ecliptic Equinox of J2000 (ECLIPJ2000) axes could not be defined" *
            " w.r.t $pname axes. Only `ICRF` or `MEME2000` are accepted as parent axes.")
        )
    end

    add_axes_inertial!(frames, axes; parent=parent, dcm=dcm)

end


"""
    add_axes_mememod!(frames, axes, parent::AbstractFrameAxes, model::IAU2006Model=iau2006b)

Add `axes` as a set of projected axes representing the Mean of Date Ecliptic Equinox to 
`frames`. 

!!! note
    Despite this class of axes has a rotation matrix that depends on time, its derivatives 
    are assumed null.

!!! warning 
    The `parent` set of axes must be the International Celestial Reference Frame (ICRF). 
    If the `parent` set of axes is not ICRF, an error is thrown.

"""
function add_axes_mememod!(frames::FrameSystem{O, T}, axes::AbstractFrameAxes,
            parent::AbstractFrameAxes) where {T, O}

    pname = axes_name(parent)
    if pname != :ICRF 
        throw(
            ErrorException("Mean Equator, Mean Equinox of date axes can only be defined " * 
            "w.r.t the International Celestial Reference Frame (ICRF)")
        )
    end

    add_axes_projected!(frames, axes, parent, Orient.orient_icrf_to_mememod)
end
