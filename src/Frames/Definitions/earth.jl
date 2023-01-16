export add_axes_itrf!

"""
    add_axes_itrf!(frames, axes, parent, model::Orient.IAU2006Model=Orient.iau2006b) 

Add `axes` as a set of axes representing the International Terrestrial Reference Frame (ITRF)
to `frames`. Only the :ICRF and :GCRF are accepted as `parent` axes. Use the `model` argument 
to specify which IAU2010 convention (IAU2006A or IAU2000B) model should be used for the 
computations. The default is set to `iau2006b`.
"""
function add_axes_itrf!(frames::FrameSystem{O, T}, 
    axes::AbstractFrameAxes, parent::AbstractFrameAxes, 
    model::Orient.IAU2006Model=Orient.iau2006b) where {T, O}

    pname = axes_name(parent)

    if !(pname in :ICRF, :GCRF) && axes_id(axes) != Orient.AXESID_ICRF
        throw(
            ErrorException("International Terrestrial Reference Frame (ITRF) axes could "*
            " not be defined w.r.t $pname axes. Only the `ICRF` or `GCRF` are accepted as"*
            " parent axes.")
        )
    end 

    add_axes_rotating!(frames, axes, parent, 
                        t->Orient.orient_itrf_to_gcrf(t, model),
                        t->Orient.orient_d_itrf_to_gcrf(t, model),
                        t->Orient.orient_dd_itrf_to_gcrf(t, model)
    )

end