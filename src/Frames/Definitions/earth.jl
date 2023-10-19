export add_axes_itrf!

""" 
    AXESID_ITRF 
   
NAIF Axes ID for the International Terrestrial Reference Frame (ITRF)
"""
const AXESID_ITRF = 3

"""
    add_axes_itrf!(frames, axes, parent, model::Orient.IAU2006Model=Orient.iau2006b) 

Add `axes` as a set of axes representing the International Terrestrial Reference Frame (ITRF)
to `frames`. Only the :ICRF and :GCRF are accepted as `parent` axes. Use the `model` argument 
to specify which IAU model model should be used for the computations. The default is set to 
`iau2006b`.

### See also 
See also [`orient_rot3_itrf_to_gcrf`](@ref)

----

    add_axes_itrf!(frames::FrameSystem{O,T}, name::Symbol, parentid::Int, 
        model::Orient.IAU2006Model=Orient.iau2006b) where {T,O}

Add a new inertial axes representing the International Terrestrial Reference Frame (ITRF)
giving a `name` and a `parent`. 
The axesid is automatically assigned as [`AXESID_ITRF`](@ref).

"""
function add_axes_itrf!(
    frames::FrameSystem{O,T},
    name::Symbol,
    parentid::Int,
    model::Orient.IAUModel=Orient.iau2006b,
) where {T,O}

    if !(parentid in (Orient.AXESID_ICRF, AXESID_GCRF))
        throw(
            ArgumentError(
                "International Terrestrial Reference Frame (ITRF) axes could " *
                "not be defined w.r.t $parentid axes. Only the `ICRF` or `GCRF` are accepted as " *
                "parent axes.",
            ),
        )
    end

    return add_axes_rotating!(
        frames,
        name,
        AXESID_ITRF,
        parentid,
        t -> adjoint(Orient.orient_rot3_itrf_to_gcrf(model, t)),
        t -> adjoint.(Orient.orient_rot6_itrf_to_gcrf(model, t)),
        t -> adjoint.(Orient.orient_rot9_itrf_to_gcrf(model, t)),
        t -> adjoint.(Orient.orient_rot12_itrf_to_gcrf(model, t)),
    )
end

function add_axes_itrf!(
    frames::FrameSystem{O,T},
    axes::AbstractFrameAxes,
    parent::AbstractFrameAxes,
    model::Orient.IAUModel=Orient.iau2006b,
) where {T,O}
    pname = axes_name(parent)

    if !(pname in (:ICRF, :GCRF)) && axes_id(parent) != Orient.AXESID_ICRF
        throw(
            ArgumentError(
                "International Terrestrial Reference Frame (ITRF) axes could " *
                "not be defined w.r.t $pname axes. Only the `ICRF` or `GCRF` are accepted as " *
                "parent axes.",
            ),
        )
    end

    return add_axes_rotating!(
        frames,
        axes,
        parent,
        t -> adjoint(Orient.orient_rot3_itrf_to_gcrf(model, t)),
        t -> adjoint.(Orient.orient_rot6_itrf_to_gcrf(model, t)),
        t -> adjoint.(Orient.orient_rot9_itrf_to_gcrf(model, t)),
        t -> adjoint.(Orient.orient_rot12_itrf_to_gcrf(model, t)),
    )
end
