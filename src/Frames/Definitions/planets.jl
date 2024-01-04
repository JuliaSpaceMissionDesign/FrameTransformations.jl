export add_axes_bcrtod!, add_axes_bci2000!

"""
    add_axes_bcrtod!(frames, axes::AbstractFrameAxes, center::AbstractFramePoint, data)

Add `axes` as a set of Body-Centered Rotating (BCR), True-of-Date (TOD) axes to the 
`frames` system. The center point (i.e., the reference body) is `center`. `data` is a dictionary 
containing a parsed `TPC` file. These axes are the equivalent of SPICE's `IAU_<BODY_NAME>` frames.

!!! warning 
    The parent axes are automatically set to the ICRF (ID = $(Orient.AXESID_ICRF)). If the 
    ICRF is not defined in `frames`, an error is thrown.

----

    add_axes_bcrtod!(frames, name::Symbol, axesid::Int, cname::Symbol, cid::Int, data)

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) type 
via the [`@axes`](@ref) macro and of an [`AbstractFramePoint`](@ref) via the [`@point`](@ref)
macro.

### See also 
See also [`add_axes_rotating!`](@ref), [`add_axes_bci2000!`](@ref) and [`Orient.AXESID_ICRF`](@ref).
"""
function add_axes_bcrtod!(
    frames::FrameSystem, axes::AbstractFrameAxes, center::AbstractFramePoint, data
)
    return add_axes_bcrtod!(
        frames, axes_name(axes), axes_id(axes), point_name(center), point_id(center), data
    )
end

# Low-level function
function add_axes_bcrtod!(
    frames::FrameSystem, 
    name::Symbol, 
    axesid::Int, 
    cname::Symbol, 
    cid::Int, 
    data::AbstractDict
)

    if !(has_axes(frames, Orient.AXESID_ICRF))
        throw(
            ErrorException(
                "Body-Centered Rotating (BCR), True-of-Date (TOD) axes can only be defined" * 
                " w.r.t. the ICRF (ID = $(Orient.AXESID_ICRF)), which is not defined in" * 
                " the current frames graph."
            )
        )
    end

    # Get nutation - precession angles 
    p = Orient.PlanetsPrecessionNutation(cid, data)

    # Build transformations 
    f1, f2, f3 = _get_bcrtod_funs(p, cname)
    
    # Insert new axes in the frame system 
    return add_axes_rotating!(frames, name, axesid, Orient.AXESID_ICRF, f1, f2, f3)

end

# Lowest-level function 
function _get_bcrtod_funs(p::Orient.PlanetsPrecessionNutation, bodyname::Symbol)
    f1, f2, f3 = Orient.orient_planets_angles(p, bodyname)

    fname = Symbol("orient_rot3_icrf_to_bcr_tod_$bodyname")
    fdname = Symbol("orient_rot6_icrf_to_bcr_tod_$bodyname")
    fddname = Symbol("orient_rot9_icrf_to_bcr_tod_$bodyname")

    return eval(
        quote
            (
                function ($fname)(sec::Number)
                    ra, dec, w = $f1(sec)
                    return angle_to_dcm(π / 2 + ra, π / 2 - dec, w, :ZXZ)
                end,

                function ($fdname)(sec::Number)
                    ra, dec, w, rad, decd, wd = $f2(sec)

                    R = angle_to_dcm(π / 2 + ra, π / 2 - dec, w, :ZXZ)
                    dR = angle_to_δdcm(
                        (π / 2 + ra, rad), (π / 2 - dec, -decd), (w, wd), :ZXZ
                    )

                    return R, dR
                end,

                function ($fddname)(sec::Number)
                    ra, dec, w, rad, decd, wd, radd, decdd, wdd = $f3(sec)

                    R = angle_to_dcm(π / 2 + ra, π / 2 - dec, w, :ZXZ)
                    dR = angle_to_δdcm(
                        (π / 2 + ra, rad), (π / 2 - dec, -decd), (w, wd), :ZXZ
                    )

                    d²R = angle_to_δ²dcm(
                        (π / 2 + ra, rad, radd),
                        (π / 2 - dec, -decd, -decdd),
                        (w, wd, wdd),
                        :ZXZ,
                    )

                    return R, dR, d²R
                end,
            )
        end,
    )
end


"""
    add_axes_bci2000!(frames, axes::AbstractFrameAxes, center, data)

Add `axes` as a set of Body-Centered Inertial (BCI) axes at J2000 relative to the `frames` 
system. The center point (i.e., the reference body) is `center` and can either be the point 
ID or an [`AbstractFramePoint`](@ref) instance. `data` is a dictionary containing a parsed 
`TPC` file. 

!!! warning 
    The parent axes are automatically set to the ICRF (ID = $(Orient.AXESID_ICRF)). If the 
    ICRF is not defined in `frames`, an error is thrown.

----

    add_axes_bci2000!(frames, name::Symbol, axesid::Int, cid::Int, data)

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) type 
via the [`@axes`](@ref) macro.
    
### See also 
See also [`add_axes_fixedoffset!`](@ref), [`add_axes_bcrtod!`](@ref) and [`Orient.AXESID_ICRF`](@ref).
"""
function add_axes_bci2000!(frames::FrameSystem, axes::AbstractFrameAxes, center, data)
    return add_axes_bci2000!(
        frames, axes_name(axes), axes_id(axes), point_alias(center), data
    )
end

# Low-level function
function add_axes_bci2000!(
    frames::FrameSystem,
    name::Symbol, 
    axesid::Int, 
    cid::Int, 
    data::AbstractDict
)

    if !(has_axes(frames, Orient.AXESID_ICRF))
        throw(
            ErrorException(
                "Body-Centered Inertial (BCI) at J2000 axes can only be defined" * 
                " w.r.t. the ICRF (ID = $(Orient.AXESID_ICRF)), which is not defined in" *
                " the current frames graph."
            )
        )
    end

    # Get nutation - precession angles 
    p = Orient.PlanetsPrecessionNutation(cid, data)

    # Build transformation
    dcm = _get_bci2000_rotation(p)

    # Insert new axes in the frame system 
    return add_axes_fixedoffset!(frames, name, axesid, Orient.AXESID_ICRF, dcm)
    
end

# Lowest-level function
function _get_bci2000_rotation(p::Orient.PlanetsPrecessionNutation)

    angles = Orient._iau_angles(p)
    a, del = angles[1][1], angles[2][1]
    Θᵢ = angles[1][5]

    @eval begin
        T = 0.0
        d = 0.0
        Θᵢ = $Θᵢ
        a2000 = $a
        d2000 = $del
        dcm = angle_to_dcm(π / 2 + a2000, π / 2 - d2000, :ZX)
    end

    return dcm
end
