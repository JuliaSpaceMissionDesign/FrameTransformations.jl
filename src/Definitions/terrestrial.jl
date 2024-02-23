export add_axes_itrf!, add_axes_mod!, add_axes_tod!, add_axes_pef!

"""
    add_axes_itrf!(frames, axes::AbstractFrameAxes, parent, model::IERSConventions.IERSModel=iers2010b) 

Add `axes` as a set of axes representing the International Terrestrial Reference Frame (ITRF)
to `frames`. Use the `model` argument to specify which IERS convention should be used for 
the computations. The default is set to `iers2010b`.

!!! warning 
    If the ID of the parent set of `axes` is neither the ICRF (ID = $(AXESID_ICRF))
    nor the GCRF (ID = $(AXESID_GCRF)), an error is thrown. 

----

    add_axes_itrf!(frames, name::Symbol, parentid::Int, conv::IERSConventions.IERSModel=iers2010b, 
        axesid::Int = AXESID_ITRF)

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) type 
via the [`@axes`](@ref) macro.

### See also 
See also [`add_axes_rotating!`](@ref).
"""
@inline function add_axes_itrf!(
    frames::FrameSystem,
    axes::AbstractFrameAxes,
    parent,
    model::IERSConventions.IERSModel=iers2010b,
)
    return add_axes_itrf!(frames, axes_name(axes), axes_alias(parent), model, axes_id(axes))
    
end

# Low-level function
function add_axes_itrf!(
    frames::FrameSystem,
    name::Symbol,
    parentid::Int,
    conv::IERSConventions.IERSModel=iers2010b,
    axesid::Int=AXESID_ITRF
)

    if !(parentid in (AXESID_ICRF, AXESID_GCRF))
        throw(
            ArgumentError(
                "International Terrestrial Reference Frame (ITRF) axes cannot be defined" *
                "with respect to axes $parentid. Only the `ICRF` (ID = $(AXESID_ICRF)) " *
                "or the `GCRF` (ID = $(AXESID_GCRF)) are accepted as parent axes.",
            ),
        )
    end

    return add_axes_rotating!(
        frames,
        name,
        axesid,
        parentid,
        t -> iers_rot3_gcrf_to_itrf(t, conv),
        t -> iers_rot6_gcrf_to_itrf(t, conv),
        t -> iers_rot9_gcrf_to_itrf(t, conv),
        t -> iers_rot12_gcrf_to_itrf(t, conv),
    )
end

"""
    add_axes_mod!(frames, axes::AbstractFrameAxes, parent, 
        model::IERSConventions.IERSModel=iers2010b)

Add `axes` as a set of projected axes representing the Mean Equator and Equinox of Date (MOD)
to `frames`. 

!!! note
    Despite this class of axes has a rotation matrix that depends on time, its derivatives 
    are assumed null.

!!! warning 
    The ID of the `parent` set of axes must be $(AXESID_ICRF) (ICRF), or $(AXESID_GCRF)
    (GCRF) otherwise an error is thrown. 

----

    add_axes_mod!(frames, name::Symbol, axesid::Int, parentid::Int, 
        model::IERSConventions.IERSModel=iers2010b)

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) type 
via the [`@axes`](@ref) macro.

### See also 
See also [`add_axes_projected!`](@ref)
"""
@inline function add_axes_mod!(frames::FrameSystem, axes::AbstractFrameAxes, parent, 
    model::IERSConventions.IERSModel=iers2010b)
    return add_axes_mod!(frames, axes_name(axes), axes_id(axes), axes_alias(parent), model)
end

# Low-level function
function add_axes_mod!(frames::FrameSystem, name::Symbol, axesid::Int, 
    parentid::Int=AXESID_ICRF, model::IERSConventions.IERSModel=iers2010b)

    if parentid != AXESID_ICRF || parentid != AXESID_GCRF
        throw(
            ArgumentError(
                "Mean of Date (MOD) axes can only be defined " *
                "w.r.t. the ICRF (ID = $(AXESID_ICRF)) or " * 
                "the GCRF (ID = $(AXESID_GCRF)).",
            )
        )
    end

    return add_axes_projected!(
        frames, name, axesid, parentid, 
        t->iers_rot3_gcrf_to_mod(t, model)
    )
end

"""
    add_axes_tod!(frames, axes::AbstractFrameAxes, parent, 
        model::IERSConventions.IERSModel=iers2010b)

Add `axes` as a set of projected axes representing the True Equator of Date (TOD)
to `frames`. 

!!! note
    Despite this class of axes has a rotation matrix that depends on time, its derivatives 
    are assumed null.

!!! warning 
    The ID of the `parent` set of axes must be $(AXESID_ICRF) (ICRF), or $(AXESID_GCRF)
    (GCRF) otherwise an error is thrown. 

----

    add_axes_tod!(frames, name::Symbol, axesid::Int, parentid::Int, 
        model::IERSConventions.IERSModel=iers2010b)

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) type 
via the [`@axes`](@ref) macro.

### See also 
See also [`add_axes_projected!`](@ref).
"""
@inline function add_axes_tod!(frames::FrameSystem, axes::AbstractFrameAxes, parent, 
    model::IERSConventions.IERSModel=iers2010b)
    return add_axes_mod!(frames, axes_name(axes), axes_id(axes), axes_alias(parent), model)
end

# Low-level function
function add_axes_tod!(frames::FrameSystem, name::Symbol, axesid::Int, 
    parentid::Int=AXESID_ICRF, model::IERSConventions.IERSModel=iers2010b)

    if parentid != AXESID_ICRF || parentid != AXESID_GCRF
        throw(
            ArgumentError(
                "True of Date (TOD) axes can only be defined " *
                "w.r.t. the ICRF (ID = $(AXESID_ICRF)) or " * 
                "the GCRF (ID = $(AXESID_GCRF)).",
            )
        )
    end

    return add_axes_projected!(
        frames, name, axesid, parentid, 
        t -> iers_rot3_gcrf_to_tod(t, model)
    )
end

"""
    add_axes_pef!(frames, axes::AbstractFrameAxes, parent,
        model::IERSConventions.IERSModel=iers2010b)

Add `axes` as a set of rotating axes representing the Pseudo Earth Fixed (PEF) axes
to `frames`. 

!!! warning 
    The ID of the `parent` set of axes must be $(AXESID_ICRF) (ICRF), 
    or $(AXESID_GCRF) (GCRF) otherwise an error is thrown. 

----

    add_axes_pef!(frames, name::Symbol, axesid::Int, parentid::Int, 
        model::IERSConventions.IERSModel=iers2010b

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) type 
via the [`@axes`](@ref) macro.

### See also 
See also [`add_axes_rotating!`](@ref).
"""
@inline function add_axes_pef!(
    frames::FrameSystem, axes::AbstractFrameAxes, parent, 
    model::IERSConventions.IERSModel=iers2010b)
    return add_axes_pef!(frames, axes_name(axes), axes_id(axes), axes_alias(parent), model)
end

# Low-level function
function add_axes_pef!(frames::FrameSystem, name::Symbol, axesid::Int, 
    parentid::Int=Orient.AXESID_ICRF, model::IERSConventions.IERSModel=iers2010b)

    if parentid != AXESID_ICRF || parentid != AXESID_GCRF
        throw(
            ArgumentError(
                "Pseudo Earth Fixed (PEF) axes can only be defined " *
                "w.r.t. the ICRF (ID = $(AXESID_ICRF)) or " * 
                "the GCRF (ID = $(AXESID_GCRF)).",
            )
        )
    end

    # TODO: PEF axes use AD for higher-order derivatives
    return add_axes_rotating!(
        frames,
        name,
        axesid,
        parentid,
        t -> iers_rot3_gcrf_to_pef(model, t)
    )
end