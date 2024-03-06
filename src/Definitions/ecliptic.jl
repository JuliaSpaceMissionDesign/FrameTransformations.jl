export add_axes_ecl2000!


"""
    add_axes_ecl2000!(frames, axes::AbstractFrameAxes, parent, 
        model::IERSConventions.IERSModel=iers1996)
    
Add `axes` as a set of inertial axes representing the Ecliptic Equinox of J2000 (ECL2000)
to `frames`. The obliquity of the ecliptic is computed using the IERS convention `model`.

The admissed `parent` set of axes are the following: 
- **ICRF**: for the International Celestial Reference Frame, with ID = $(AXESID_ICRF)
- **GCRF**: for the Geocentric Celestial Reference Frame, with ID = $(AXESID_GCRF)
- **EME2000**: the Mean Earth/Moon Ephemeris of J2000, with ID = $(AXESID_EME2000)

!!! note 
    To retrieve the same orientation of the ECLIPJ2000 axes avaialble in the SPICE 
    Toolkit, the `iers1996` model must be used.

----

    add_axes_ecl2000!(frames, name::Symbol, parentid::Int, 
        model::IERSConventions.IERSModel=iers1996, axesid::Int = AXESID_ECL2000)

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) type 
via the [`@axes`](@ref) macro.

### Examples
```julia-repl 
julia> @axes ICRF 1 InternationalCelestialReferenceFrame

julia> @axes ECL2000 17 EclipticEquinoxJ2000 

julia> FRAMES = FrameSystem{1, Float64}();

julia> add_axes_inertial!(FRAMES, ICRF)

julia> add_axes_ecl2000!(FRAMES, ECL2000, ICRF)
``` 

### See also 
See also [`add_axes_inertial!`](@ref) and [`DCM_ICRF_TO_ECL2000`](@ref)
"""
@inline function add_axes_ecl2000!(
    frames::FrameSystem,
    axes::AbstractFrameAxes,
    parent,
    model::IERSConventions.IERSModel=iers1996,
)
    return add_axes_ecl2000!(
        frames, axes_name(axes), axes_alias(parent), model, axes_id(axes)
    )

end

# Low-level function
function add_axes_ecl2000!(
    frames::FrameSystem, 
    name::Symbol, 
    parentid::Int, 
    model::IERSConventions.IERSModel=iers1996,
    axesid::Int = AXESID_ECL2000
)

    DCM_EME2000_TO_ECLJ2000 = angle_to_dcm(iers_obliquity(model, 0.0), :X)

    if parentid == AXESID_ICRF || parentid == AXESID_GCRF
        dcm = DCM_EME2000_TO_ECLJ2000_ * DCM_ICRF_TO_EME2000
    elseif parentid == AXESID_EME2000
        # Compute the J2000 to ECLIPJ2000 rotationa ccording to the desired IAU model
        dcm = DCM_EME2000_TO_ECLJ2000
    else 
        throw(
            ArgumentError(
                "Ecliptic Equinox of J2000 (ECL2000) axes cannot be defined" *
                " w.r.t. $parentid axes. Only `ICRF` (ID = $(AXESID_ICRF)) or" * 
                " `EME2000` (ID = $(AXESID_EME2000)) are accepted as parent axes.",
            ),
        )
    end

    if axesid != AXESID_ECL2000
        @warn "$name is aliasing an ID that is not the standard ECL2000 ID" *
              " ($(AXESID_ECL2000))."
    end

    return add_axes_inertial!(
        frames, name, axesid; parentid = parentid, dcm = dcm)

end
