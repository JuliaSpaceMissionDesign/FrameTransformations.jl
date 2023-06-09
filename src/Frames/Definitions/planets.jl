export add_axes_bcrtod!, add_axes_bci2000!

"""
    add_axes_bcrtod!(frames::FrameSystem{O, T}, 
        data::AbstractDict, center::AbstractFramePoint, 
        axes::AbstractFrameAxes, parent::AbstractFrameAxes) where {T, O}

Insert a Body-Centered Rotating (BCR), True-of-Date (TOD) axes to the `frames` system.

### Input/s

- `frames` - The frame system to which the new frame will be added.
- `data` - A dictionary containing a parsed `TPC` file. The dictionary as has keys the NAIFId 
    of the bodies contained in the file, and as values a dictionary containing the keywords 
    and the values of the file's data.
- `center`: The center point of the new axes.
- `axes`: The new axes to be added to the frame system.
- `parent`: The parent axes of the new ones.

!!! note 
    The axes constructed here corresponds to the SPICE `IAU_<BODY_NAME>` frames. 

!!! warning 
    The `parent` set of axes must be the International Celestial Reference Frame (ICRF). 
    If the `parent` set of axes is not ICRF, an error is thrown.
"""
function add_axes_bcrtod!(
    frames::FrameSystem{O,T},
    data::AbstractDict,
    center::AbstractFramePoint,
    axes::AbstractFrameAxes,
    parent::AbstractFrameAxes,
) where {T,O}
    pname = axes_name(parent)
    if pname != :ICRF
        throw(
            ErrorException(
                "Body-Centered Rotating (TOD) axes could be defined only" *
                " w.r.t the International Celestial Reference Frame (ICRF)",
            ),
        )
    end

    return _axes_bcrtod!(
        frames, data, point_id(center), String(point_name(center)), axes, parent
    )
end

function _axes_bcrtod!(
    frames::FrameSystem{O,T},
    data::AbstractDict,
    cid::Integer,
    cname::String,
    axes::AbstractFrameAxes,
    parent::AbstractFrameAxes,
) where {T,O}

    # Get nutation - precession angles 
    p = Orient.PlanetsPrecessionNutation(cid, data)

    # Build transformations 
    f1, f2, f3 = _axes_bcrtod(p, cname)

    # Insert new axes in the frame system 
    return add_axes_rotating!(frames, axes, parent, f1, f2, f3)
end

function _axes_bcrtod(p::Orient.PlanetsPrecessionNutation, bodyname)
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
    add_axes_bci2000!(frames::FrameSystem{O, T}, 
        data::AbstractDict, center::AbstractFramePoint, 
        axes::AbstractFrameAxes, parent::AbstractFrameAxes) where {T, O}

Insert Body-Centered Inertial (BCI) axes at J2000 relative to the body `center` to the 
`frames` system.

### Input/s

- `frames` - The frame system to which the new frame will be added.
- `data` - A dictionary containing a parsed `TPC` file. The dictionary as has keys the NAIFId 
    of the bodies contained in the file, and as values a dictionary containing the keywords 
    and the values of the file's data.
- `center`: The center point of the new axes.
- `axes`: The new axes to be added to the frame system.
- `parent`: The parent axes of the new ones.

!!! warning 
    The `parent` set of axes must be the International Celestial Reference Frame (ICRF). 
    If the `parent` set of axes is not ICRF, an error is thrown.
"""
function add_axes_bci2000!(
    frames::FrameSystem{O,T},
    data::AbstractDict,
    center::AbstractFramePoint,
    axes::AbstractFrameAxes,
    parent::AbstractFrameAxes,
) where {T,O}
    pname = axes_name(parent)
    if pname != :ICRF
        throw(
            ErrorException(
                "Body-Centered Inertial axes could be defined only" *
                " w.r.t the International Celestial Reference Frame (ICRF)",
            ),
        )
    end

    return _axes_bci2000!(frames, data, point_id(center), axes, parent)
end

function _axes_bci2000!(
    frames::FrameSystem{O,T},
    data::AbstractDict,
    cid::Integer,
    axes::AbstractFrameAxes,
    parent::AbstractFrameAxes,
) where {T,O}

    # Get nutation - precession angles 
    p = Orient.PlanetsPrecessionNutation(cid, data)

    # Build transformation
    # In this case the transformation is a static rotation matrix 
    dcm = _axes_bci2000(p)

    # Insert new axes in the frame system 
    return add_axes_fixedoffset!(frames, axes, parent, dcm)
end

function _axes_bci2000(p::Orient.PlanetsPrecessionNutation)
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
