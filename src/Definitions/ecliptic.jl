
function add_axes_ecl2000!(
    frames::FrameSystem, name::Symbol, parentid::Int=AXESID_ICRF, id::Int=AXESID_ECL2000,
    model::IERSModel=iers2010a,
)

    # Compute the J2000 to ECLIPJ2000 rotation according to the desired IERS model
    DCM_EME2000_TO_ECLJ2000 = angle_to_dcm(iers_obliquity(model, 0), :X)

    if parentid == AXESID_ICRF || parentid == AXESID_GCRF
        dcm = DCM_EME2000_TO_ECLJ2000 * DCM_ICRF_TO_EME2000
    elseif parentid == AXESID_EME2000
        dcm = DCM_EME2000_TO_ECLJ2000
    else 
        throw(
            ArgumentError(
                "Ecliptic Equinox of J2000 (ECL2000) axes cannot be defined" *
                " w.r.t. $parentid axes. Only the ICRF (ID = $(AXESID_ICRF)), the GCRF" * 
                " (ID = $(AXESID_GCRF)) or the EME2000 (ID = $(AXESID_EME2000)) are" * 
                " accepted as parent axes.",
            ),
        )
    end

    if id != AXESID_ECL2000
        @warn "$name is aliasing an ID that is not the standard ECL2000 ID" *
              " ($(AXESID_ECL2000))."
    end

    return add_axes_inertial!(frames, name, id, parentid, dcm)
end