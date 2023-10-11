export add_axes_icrf!, add_axes_gcrf!

"""
    add_axes_icrf!(frames::FrameSystem)

Add the International Celestial Reference Frame (ICRF) to the frames graph. 

!!! warning 
    Shall be the root axes.
"""
function add_axes_icrf!(frames::FrameSystem)
    return add_axes_inertial!(frames, :ICRF, Orient.AXESID_ICRF)
end


"""
    add_axes_gcrf!(frames::FrameSystem; root=false)

Add the Geocentric Celestial Reference Frame (GCRF) to the frames graph. Could be used as 
`root` axes or not. Could be child only of the International Celestial Reference Frame (ICRF) 
axes.
"""
function add_axes_gcrf!(frames::FrameSystem; root=false)
    if root 
        return add_axes_inertial!(frames, :GCRF, Orient.AXESID_GCRF)
    else 
        return add_axes_fixedoffset!(frames, :GCRF, Orient.AXESID_GCRF, Orient.AXESID_ICRF, DCM(1.0I)) 
    end
end

