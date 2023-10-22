export add_axes_icrf!, add_axes_gcrf!

"""
    add_axes_icrf!(frames::FrameSystem)

Add the International Celestial Reference Frame (ICRF) as the root axes of the frames graph.
The axes are automatically named `ICRF` and assigned the $(Orient.AXESID_ICRF) ID. 

### See also 
See also [`add_axes_inertial!`](@ref), [`add_axes_gcrf!`](@ref) and [`Orient.AXESID_ICRF`](@ref).
"""
@inline function add_axes_icrf!(frames::FrameSystem)

    if !isempty(frames_axes(frames))
        throw(ArgumentError("The ICRF can only be defined as a set of root axes."))
    end

    return add_axes_inertial!(frames, :ICRF, Orient.AXESID_ICRF)

end


"""
    add_axes_gcrf!(frames::FrameSystem)

Add the Geocentric Celestial Reference Frame (GCRF) to the frames graph. The axes are 
automatically named `GCRF` and assigned the $((Orient.AXESID_GCRF)) ID. These axes can only 
be defined as a set of root axes or as child of the ICRF (ID = $(Orient.AXESID_ICRF)).

### See also 
See also [`add_axes_inertial!`](@ref), [`add_axes_icrf!`](@ref) and [`Orient.AXESID_GCRF`](@ref).
"""
function add_axes_gcrf!(frames::FrameSystem)

    if has_axes(frames, Orient.AXESID_ICRF)
        # Add the GCRF as a child of the ICRF with an identity rotation 
        return add_axes_fixedoffset!(
            frames, :GCRF, Orient.AXESID_GCRF, Orient.AXESID_ICRF, DCM(1.0I)
        )

    elseif isempty(frames_axes(frames))
        # Add the GCRF as a root set of axes
        return add_axes_inertial!(frames, :GCRF, Orient.AXESID_GCRF)
        
    else 
        throw(
            ArgumentError(
                "The GCRF can only be defined with respect to the ICRF (ID =" * 
                " $(Orient.AXESID_ICRF)) or as a set of root axes."
            )
        )
    end
end

