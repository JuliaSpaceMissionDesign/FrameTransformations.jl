export DCM_ICRF_TO_EME2000, DCM_EME2000_TO_ECL2000, DCM_ICRF_TO_ECL2000,
       add_axes_icrf!, add_axes_gcrf!, add_axes_eme2000!

"""
    DCM_ICRF_TO_EME2000

DCM for the rotation from the International Celestial Reference Frame (`ICRF`) and the 
Mean Equator and Equinox of J2000.0 (`EME2000`). This corresponds to the `J2000` frame in 
the SPICE toolkit.

!!! note 
    The frame bias is here computed using the IAU 2006 Precession model, similarly to ESA's 
    GODOT. Some other software libraries, such as Orekit, use the frame bias of the IAU 2000 
    precession model. The two definitions differ of about 1 arcsecond.

    Moreover, according to [Hilton](https://www.aanda.org/articles/aa/pdf/2004/02/aa3851.pdf) 
    there are multiple possibilities to define the proper rotation between the ICRS and 
    the EME2000. The transformation implemented here correspond to Eq. 6 using the parameters 
    in Table 3, line 1 (RIERS).

### References
- Hilton, James L., and Catherine Y. Hohenkerk. -- Rotation matrix from the mean 
    dynamical equator and equinox at J2000. 0 to the ICRS. -- Astronomy & Astrophysics 
    513.2 (2004): 765-770. DOI: [10.1051/0004-6361:20031552](https://www.aanda.org/articles/aa/pdf/2004/02/aa3851.pdf)
- [SOFA docs](https://www.iausofa.org/2021_0512_C/sofa/sofa_pn_c.pdf)
"""
const DCM_ICRF_TO_EME2000 = IERSConventions.iers_pb(iers2010a, 0.0)

"""
    DCM_EME2000_TO_ECL2000

DCM for the rotation from the Mean Equator and Equinox of J2000 (`EME2000`) to the 
Mean Ecliptic Equinox. This corresponds to the transformation `J2000 -> ECLIPJ2000` 
in the SPICE toolkit, and uses the mean obliquity of the ecliptic from the IAU 1976 theory.
"""
const DCM_EME2000_TO_ECL2000 = angle_to_dcm(iers_obliquity(iers1996, 0.0), :X)

"""
    DCM_ICRF_TO_ECL2000

DCM for the rotation from the International Celestial Reference Frame (`ICRF`) to the 
Mean Ecliptic Equinox of J2000 (`ECL2000`).
"""
const DCM_ICRF_TO_ECL2000 = DCM_ICRF_TO_EME2000 * DCM_EME2000_TO_ECL2000


"""
    add_axes_icrf!(frames::FrameSystem)

Add the International Celestial Reference Frame (ICRF) as the root axes of the frames graph.
The axes are automatically named `ICRF` and assigned the $(AXESID_ICRF) ID. 

### See also 
See also [`add_axes_inertial!`](@ref), [`add_axes_gcrf!`](@ref) and [`AXESID_ICRF`](@ref).
"""
@inline function add_axes_icrf!(frames::FrameSystem)
    if !isempty(frames_axes(frames))
        throw(ArgumentError("The ICRF can only be defined as a set of root axes."))
    end

    return add_axes_inertial!(frames, :ICRF, AXESID_ICRF)
end

"""
    add_axes_gcrf!(frames::FrameSystem)

Add the Geocentric Celestial Reference Frame (GCRF) to the frames graph. The axes are 
automatically named `GCRF` and assigned the $((AXESID_GCRF)) ID. These axes can only 
be defined as a set of root axes or as child of the ICRF (ID = $(AXESID_ICRF)).

### See also 
See also [`add_axes_inertial!`](@ref), [`add_axes_icrf!`](@ref) and [`AXESID_GCRF`](@ref).
"""
function add_axes_gcrf!(frames::FrameSystem)
    if has_axes(frames, AXESID_ICRF)
        # Add the GCRF as a child of the ICRF with an identity rotation 
        return add_axes_fixedoffset!(
            frames, :GCRF, AXESID_GCRF, AXESID_ICRF, DCM(1.0I)
        )

    elseif isempty(frames_axes(frames))
        # Add the GCRF as a root set of axes
        return add_axes_inertial!(frames, :GCRF, AXESID_GCRF)
        
    else 
        throw(
            ArgumentError(
                "The GCRF can only be defined with respect to the ICRF (ID =" * 
                " $(AXESID_ICRF)) or as a set of root axes."
            )
        )
    end
end


"""
    add_axes_eme2000!(frames, axes::AbstractFrameAxes, parent)
    
Add `axes` as a set of inertial axes representing the Mean Equator Mean Equinox of J2000 
to `frames`. 

!!! warning 
    The the axes ID of the parent set of axes must be $(AXESID_ICRF) (ICRF), 
    $(AXESID_GCRF) (GCRF) or $(AXESID_ECL2000) (EME2000) otherwise and error is thrown.

----

    add_axes_eme2000!(frames, name::Symbol, parentid::Int, axesid::Int = AXESID_EME2000)

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) type 
via the [`@axes`](@ref) macro.

### Examples 
```julia-repl 
julia> @axes ICRF 1 InternationalCelestialReferenceFrame

julia> @axes EME2000 22 MeanEquatorMeanEquinoxJ2000 

julia> FRAMES = FrameSystem{1, Float64}();

julia> add_axes_inertial!(FRAMES, ICRF)

julia> add_axes_eme2000!(FRAMES, EME2000, ICRF)
```

### See also 
See also [`add_axes_inertial!`](@ref) and [`DCM_ICRF_TO_EME2000`](@ref)
"""
@inline function add_axes_eme2000!(frames::FrameSystem, axes::AbstractFrameAxes, parent)
    return add_axes_eme2000!(frames, axes_name(axes), axes_alias(parent), axes_id(axes))
end

# Low-level function
function add_axes_eme2000!(
    frames::FrameSystem, name::Symbol, parentid::Int=AXESID_ICRF, 
    axesid::Int=AXESID_EME2000
)
    if parentid == AXESID_ICRF || parentid == AXESID_GCRF
        dcm = DCM_ICRF_TO_EME2000
    elseif parentid == AXESID_ECL2000
        dcm = DCM_EME2000_TO_ECL2000'
    else 
        throw(
            ArgumentError(
                "Mean Equator, Mean Equinox of J2000 (EME2000) axes can only be defined " *
                "w.r.t. the ICRF (ID = $(AXESID_ICRF)), the GCRF (ID = $(AXESID_GCRF)) " * 
                "or the ECL2000 (ID = $(AXESID_ECL2000)).",
            ),
        )
    end

    if axesid != AXESID_EME2000
        @warn "$name is aliasing an ID that is not the standard EME2000 ID" *
              " ($(AXESID_EME2000))."
    end

    return add_axes_inertial!(frames, name, axesid; parentid = parentid, dcm = dcm)
end
