
"""
    add_axes_itrf!(fr::FrameSystem, name::Symbol, parentid::Int=AXESID_ICRF, id::Int=AXESID_ITRF,
        model::IERSModel=iers2010b) 

Add International Terrestrial Reference Frame (ITRF) axes to `fr`. Use the `model` 
argument to specify which IERS convention should be used for the computations.

!!! warning 
    If the ID of the parent set of `axes` is neither the ICRF (ID = $(AXESID_ICRF))
    nor the GCRF (ID = $(AXESID_GCRF)), an error is thrown. 
"""
function add_axes_itrf!(
    fr::FrameSystem, name::Symbol, parentid::Int=AXESID_ICRF, id::Int=AXESID_ITRF,
    model::IERSModel=iers2010b
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
        fr, name, id, parentid,
        t -> iers_rot3_gcrf_to_itrf(t, model), t -> iers_rot6_gcrf_to_itrf(t, model),
        t -> iers_rot9_gcrf_to_itrf(t, model), t -> iers_rot12_gcrf_to_itrf(t, model)
    )
end

# ------------------------------------------------------------------------------------------
# Equinox-based transformations from ICRF/GCRF
# ------------------------------------------------------------------------------------------

"""
    add_axes_mod!(fr, name::Symbol, axesid::Int, parentid::Int, 
        model::IERSConventions.IERSModel=iers2010b)

Add Mean Equator and Equinox of Date (MOD) axes to `fr`. Use the `model` argument to 
specify which IERS convention should be used for the computations.

!!! note 
    The Mean-of-Date axes are obtained by applying the frame bias and precession matrix. 
    For this reason, if the IERS 1996 conventions are used, the rotation is 
    actually computed starting from the EME2000 rather than the GCRF.  

!!! note
    Despite this class of axes has a rotation matrix that depends on time, its derivatives 
    are assumed null.

!!! warning 
    The ID of the `parent` set of axes must be $(AXESID_ICRF) (ICRF) or $(AXESID_GCRF)
    (GCRF) otherwise an error is thrown. 
"""
function add_axes_mod!(fr::FrameSystem, name::Symbol, id::Int, parentid::Int=AXESID_ICRF, 
    model::IERSModel=iers2010b
)

    if !(parentid in (AXESID_ICRF, AXESID_GCRF))
        throw(
            ArgumentError(
                "Mean-of-Date (MOD) axes cannot be defined with respect to axes $parentid."*
                " Only the ICRF (ID = $(AXESID_ICRF)) or the GCRF"*
                " (ID = $(AXESID_GCRF)) are accepted as parent axes.",
            ),
        )
    end

    return add_axes_projected!(
        fr, name, id, parentid, t->iers_rot3_gcrf_to_mod(t, model)
    )
end

"""
    add_axes_tod!(fr::FrameSystem, name::Symbol, id::Int, parentid::Int=AXESID_ICRF, 
        model::IERSModel=iers2010b)

Add True Equator of Date (TOD) axes to `fr`. Use the `model` argument to specify which 
IERS convention should be used for the computations.

!!! note 
    The True-of-Date axes are obtained by applying the frame bias, precession and 
    nutation matrix. For this reason, if the IERS 1996 conventions are used, the 
    rotation is actually computed starting from the EME2000 rather than the GCRF.  

!!! note
    Despite this class of axes has a rotation matrix that depends on time, its derivatives 
    are assumed null.

!!! warning 
    The ID of the `parent` set of axes must be $(AXESID_ICRF) (ICRF) or $(AXESID_GCRF)
    (GCRF) otherwise an error is thrown. 
"""
function add_axes_tod!(fr::FrameSystem, name::Symbol, id::Int, parentid::Int=AXESID_ICRF, 
    model::IERSModel=iers2010b
)

    if !(parentid in (AXESID_ICRF, AXESID_GCRF))
        throw(
            ArgumentError(
                "True-of-Date (TOD) axes cannot be defined with respect to axes $parentid."*
                " Only the ICRF (ID = $(AXESID_ICRF)) or the GCRF"*
                " (ID = $(AXESID_GCRF)) are accepted as parent axes.",
            ),
        )
    end

    return add_axes_projected!(
        fr, name, id, parentid, t -> iers_rot3_gcrf_to_tod(t, model)
    )
end

"""
    add_axes_gtod!(fr::FrameSystem, name::Symbol, id::Int, parentid::Int=AXESID_ICRF, 
        model::IERSModel=iers2010b)

Add Greenwich True-of-Date (GTOD) axes to `fr`. Use the `model` argument to specify 
which IERS convention should be used for the computations.

!!! note
    Despite this class of axes has a rotation matrix that depends on time, its derivatives 
    are assumed null.

!!! warning 
    The ID of the `parent` set of axes must be $(AXESID_ICRF) (ICRF) or $(AXESID_GCRF)
    (GCRF) otherwise an error is thrown. 
"""
function add_axes_gtod!(
    fr::FrameSystem, name::Symbol, id::Int, parentid::Int=AXESID_ICRF, 
    model::IERSModel=iers2010b
)

    if !(parentid in (AXESID_ICRF, AXESID_GCRF))
        throw(
            ArgumentError(
                "Greenwich True-of-Date (GTOD) axes cannot be defined with respect to axes"*
                " $parentid. Only the ICRF (ID = $(AXESID_ICRF)) or the GCRF"*
                " (ID = $(AXESID_GCRF)) are accepted as parent axes.",
            ),
        )
    end

    # TODO: check if inertial is right
    return add_axes_projected!(
        fr, name, id, parentid, t -> iers_rot3_gcrf_to_gtod(t, model)
    )
end

"""
    add_axes_pef!(fr::FrameSystem, name::Symbol, id::Int, parentid::Int=AXESID_ICRF, 
        model::IERSModel=iers2010b)

Add Pseudo-Earth Fixed (PEF) axes to `fr`. Use the `model` argument to specify which 
IERS convention should be used for the computations.

!!! warning 
    The ID of the `parent` set of axes must be $(AXESID_ICRF) (ICRF) or $(AXESID_GCRF)
    (GCRF) otherwise an error is thrown. 
"""
function add_axes_pef!(
    fr::FrameSystem, name::Symbol, id::Int, parentid::Int=AXESID_ICRF, 
    model::IERSModel=iers2010b
)

    if !(parentid in (AXESID_ICRF, AXESID_GCRF))
        throw(
            ArgumentError(
                "Pseudo-Earth Fixed (PEF) axes cannot be defined with respect to"*
                " axes $parentid. Only the ICRF (ID = $(AXESID_ICRF)) or the GCRF"*
                " (ID = $(AXESID_GCRF)) are accepted as parent axes.",
            ),
        )
    end

    # TODO: PEF axes use AD for higher-order derivatives
    return add_axes_rotating!(
        fr, name, id, parentid, t -> iers_rot3_gcrf_to_pef(model, t)
    )
end

# ------------------------------------------------------------------------------------------
# CIO-based transformations from ICRF/GCRF 
# ------------------------------------------------------------------------------------------

"""
    add_axes_cirf!(fr::FrameSystem, name::Symbol, id::Int, parentid::Int=AXESID_ICRF, 
        model::IERSModel=iers2010b)

Add Celestial Intermediate Reference Frame (CIRF) axes to `fr`. Use the `model` argument 
to specify which IERS convention should be used for the computations.

!!! note
    Despite this class of axes has a rotation matrix that depends on time, its derivatives 
    are assumed null.

!!! warning 
    The ID of the `parent` set of axes must be $(AXESID_ICRF) (ICRF) or $(AXESID_GCRF)
    (GCRF) otherwise an error is thrown. 
"""
function add_axes_cirf!(fr::FrameSystem, name::Symbol, id::Int, parentid::Int=AXESID_ICRF, 
    model::IERSModel=iers2010b
)

    if !(parentid in (AXESID_ICRF, AXESID_GCRF))
        throw(
            ArgumentError(
                "The Celestial Intermediate Reference Frame (CIRF) axes cannot be defined"*
                " with respect to axes $parentid. Only the ICRF (ID = $(AXESID_ICRF)) or"*
                " the GCRF (ID = $(AXESID_GCRF)) are accepted as parent axes.",
            ),
        )
    end

    # TODO: check if inertial is right
    return add_axes_projected!(
        fr, name, id, parentid, t -> iers_rot3_gcrf_to_cirf(t, model)
    )
end

"""
    add_axes_tirf!(fr::FrameSystem, name::Symbol, id::Int, parentid::Int=AXESID_ICRF, 
        model::IERSModel=iers2010b)

Add Terrestrial Intermediate Reference Frame (TIRF) axes to `fr`. Use the `model` argument 
to specify which IERS convention should be used for the computations.

!!! note
    Despite this class of axes has a rotation matrix that depends on time, its derivatives 
    are assumed null.

!!! warning 
    The ID of the `parent` set of axes must be $(AXESID_ICRF) (ICRF) or $(AXESID_GCRF)
    (GCRF) otherwise an error is thrown. 
"""
function add_axes_tirf!(fr::FrameSystem, name::Symbol, id::Int, parentid::Int=AXESID_ICRF, 
    model::IERSModel=iers2010b
)

    if !(parentid in (AXESID_ICRF, AXESID_GCRF))
        throw(
            ArgumentError(
                "The Terrestrial Intermediate Reference Frame (TIRF) axes cannot be defined"*
                " with respect to axes $parentid. Only the ICRF (ID = $(AXESID_ICRF)) or"*
                " the GCRF (ID = $(AXESID_GCRF)) are accepted as parent axes.",
            ),
        )
    end

    # TODO: check if inertial is right
    return add_axes_projected!(
        fr, name, id, parentid, t -> iers_rot3_gcrf_to_tirf(t, model)
    )
end
