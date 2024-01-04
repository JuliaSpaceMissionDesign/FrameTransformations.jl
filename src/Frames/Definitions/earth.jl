export add_axes_itrf!

"""
    add_axes_itrf!(frames, axes::AbstractFrameAxes, parent, model::IAU2006Model=iau2006b) 

Add `axes` as a set of axes representing the International Terrestrial Reference Frame (ITRF)
to `frames`. Use the `model` argument to specify which IAU model model should be used for 
the computations. The default is set to `iau2006b`.

!!! warning 
    If the ID of the parent set of `axes` is neither the ICRF (ID = $(Orient.AXESID_ICRF))
    nor the GCRF (ID = $(Orient.AXESID_GCRF)), an error is thrown. 

----

    add_axes_itrf!(frames, name::Symbol, parentid::Int, model::IAU2006Model=iau2006b, 
        axesid::Int = Orient.AXESID_ITRF)

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) type 
via the [`@axes`](@ref) macro.

### See also 
See also [`add_axes_rotating!`](@ref) and [`orient_rot3_itrf_to_gcrf`](@ref).
"""
@inline function add_axes_itrf!(
    frames::FrameSystem,
    axes::AbstractFrameAxes,
    parent,
    model::Orient.IAUModel=Orient.iau2006b,
)
    return add_axes_itrf!(frames, axes_name(axes), axes_alias(parent), model, axes_id(axes))
    
end

# Low-level function
function add_axes_itrf!(
    frames::FrameSystem,
    name::Symbol,
    parentid::Int,
    model::Orient.IAUModel=Orient.iau2006b,
    axesid::Int=Orient.AXESID_ITRF
)

    if !(parentid in (Orient.AXESID_ICRF, Orient.AXESID_GCRF))
        throw(
            ArgumentError(
                "International Terrestrial Reference Frame (ITRF) axes cannot be defined" *
                "with respect to axes $parentid. Only the `ICRF` (ID = $(Orient.AXESID_ICRF)) " *
                "or the `GCRF` (ID = $(Orient.AXESID_GCRF)) are accepted as parent axes.",
            ),
        )
    end

    return add_axes_rotating!(
        frames,
        name,
        axesid,
        parentid,
        t -> adjoint(Orient.orient_rot3_itrf_to_gcrf(model, t)),
        t -> adjoint.(Orient.orient_rot6_itrf_to_gcrf(model, t)),
        t -> adjoint.(Orient.orient_rot9_itrf_to_gcrf(model, t)),
        t -> adjoint.(Orient.orient_rot12_itrf_to_gcrf(model, t)),
    )
end

