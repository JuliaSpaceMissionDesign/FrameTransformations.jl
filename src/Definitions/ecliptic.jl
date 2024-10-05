"""
    add_axes_ecl2000!(frames, name::Symbol=:ECL2000, parentid::Int, id::Int = AXESID_ECL2000; 
        model::IERSModel=iers1996)
    
Add Ecliptic Equinox of J2000 (ECL2000) axes to `frames`. Custom `id`, `name` and `parentid`
can be assigned by the user. The admissible parent axes are the following: 
- *ICRF*: for the International Celestial Reference Frame, with ID = $(AXESID_ICRF)
- *GCRF*: for the Geocentric Celestial Reference Frame, with ID = $(AXESID_GCRF)
- *EME2000*: the Mean Earth/Moon Ephemeris of J2000, with ID = $(AXESID_EME2000)

!!! note 
    The obliquity of the ecliptic is computed using the IERS convention `model`.
    To retrieve the same orientation of the ECLIPJ2000 axes avaialble in the SPICE 
    Toolkit, the `iers1996` model must be used.
"""
function add_axes_ecl2000!(
    frames::FrameSystem, name::Symbol=:ECL2000, parentid::Int=AXESID_ICRF, id::Int=AXESID_ECL2000;
    model::IERSModel=iers1996,
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

    return add_axes_fixedoffset!(frames, name, id, parentid, dcm)
end