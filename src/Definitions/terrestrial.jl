
function add_axes_itrf!(
    frames::FrameSystem, name::Symbol, parentid::Int=AXESID_ICRF, id::Int=AXESID_ITRF,
    model::IERSModel=iers2010b,
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
        frames, name, id, parentid,
        t -> iers_rot3_gcrf_to_itrf(t, model), t -> iers_rot6_gcrf_to_itrf(t, model),
        t -> iers_rot9_gcrf_to_itrf(t, model), t -> iers_rot12_gcrf_to_itrf(t, model)
    )
end


# Equinox-based transformations from ICRF/GCRF
# ==============================================

function add_axes_mod!(frames::FrameSystem, name::Symbol, id::Int, parentid::Int=AXESID_ICRF, 
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

    return add_axes_inertial!(
        frames, name, id, parentid, t->iers_rot3_gcrf_to_mod(t, model)
    )
end

function add_axes_tod!(frames::FrameSystem, name::Symbol, id::Int, parentid::Int=AXESID_ICRF, 
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

    return add_axes_inertial!(
        frames, name, id, parentid, t -> iers_rot3_gcrf_to_tod(t, model)
    )
end


function add_axes_gtod!(
    frames::FrameSystem, name::Symbol, id::Int, parentid::Int=AXESID_ICRF, 
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
    return add_axes_inertial!(
        frames, name, id, parentid, t -> iers_rot3_gcrf_to_gtod(t, model)
    )
end

function add_axes_pef!(
    frames::FrameSystem, name::Symbol, id::Int, parentid::Int=AXESID_ICRF, 
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
        frames, name, id, parentid, t -> iers_rot3_gcrf_to_pef(model, t)
    )
end

# CIO-based transformations from ICRF/GCRF 
# ==============================================

function add_axes_cirf!(frames::FrameSystem, name::Symbol, id::Int, parentid::Int=AXESID_ICRF, 
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
    return add_axes_inertial!(
        frames, name, id, parentid, t -> iers_rot3_gcrf_to_cirf(t, model)
    )
end

function add_axes_tirf!(frames::FrameSystem, name::Symbol, id::Int, parentid::Int=AXESID_ICRF, 
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
    return add_axes_inertial!(
        frames, name, id, parentid, t -> iers_rot3_gcrf_to_tirf(t, model)
    )
end

