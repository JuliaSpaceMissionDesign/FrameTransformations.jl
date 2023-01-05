export add_axes_bcrtod!, add_axes_bci2000!

function _orient_bcrtod(p::Orient.PlanetsPrecessionNutation, bodyname)
    f1, f2, f3 = Orient.orient_planets_angles(p, bodyname)

    fname = Symbol("orient_axes_icrf_to_bcr_tod_$bodyname")
    fdname = Symbol("orient_axes_d_icrf_to_bcr_tod_$bodyname")
    fddname = Symbol("orient_axes_dd_icrf_to_bcr_tod_$bodyname")

    return eval(
        quote
            (
                function ($fname)(sec::Number)
                    ra, dec, w = $f1(sec)
                    return angle_to_dcm(π/2 + ra, π/2 - dec, w, :ZXZ)
                end, 
                function ($fdname)(sec::Number)
                    ra, dec, w, rad, decd, wd = $f2(sec)
                    ω = SA[rad/CENTURY2SEC, decd/CENTURY2SEC, wd/DAY2SEC]
                    R = angle_to_dcm(π/2 + ra, π/2 - dec, w, :ZXZ)
                    return R, DCM(ddcm(R, ω))
                end,
                function ($fddname)(sec::Number)
                    # TODO: test derivative
                    ra, dec, w, rad, decd, wd, radd, decdd, wdd = $f3(sec)
                    ω = SA[rad/CENTURY2SEC, decd/CENTURY2SEC, wd/DAY2SEC]
                    R = angle_to_dcm(π/2 + ra, π/2 - dec, w, :ZXZ)
                    ωd = SA[
                        radd/CENTURY2SEC/CENTURY2SEC, 
                        decdd/CENTURY2SEC/CENTURY2SEC, 
                        wdd/DAY2SEC/DAY2SEC
                    ]
                    Rd = DCM(ddcm(R, ω))
                    Rdd = DCM(ddcm(R, ωd) + ddcm(Rd, ω))
                    return R, Rd, Rdd
                end
            )
        end
    )
end

function _orient_bci2000(p::Orient.PlanetsPrecessionNutation)
    Θᵢ = Orient._compute_thetas(:T, p.ra.Θ₁, p.ra.Θ₂)
    angles = Orient._iau_angles(p)
    a, del = angles[1][1], angles[2][1]

    @eval begin
        T = 0.0
        d = 0.0 
        Θᵢ = $Θᵢ
        a2000 = $a
        d2000 = $del
        dcm = angle_to_dcm(π/2 + a2000, π/2 - d2000, :ZX)
    end

    return dcm
end


function add_axes_bcrtod!(frames::FrameSystem{O, T}, 
    data::AbstractDict, center::AbstractFramePoint, 
    axes::AbstractFrameAxes, parent::AbstractFrameAxes) where {T, O}

    # Get nutation - precession angles 
    p = Orient.PlanetsPrecessionNutation(point_id(center), data)

    # Build transformations 
    f1, f2, f3 = _orient_bcrtod(p, String(point_name(center)))

    # Insert new axes in the frame system 
    add_axes_rotating!(frames, axes, parent, f1, f2)
end


function add_axes_bci2000!(frames::FrameSystem{O, T}, 
    data::AbstractDict, center::AbstractFramePoint, 
    axes::AbstractFrameAxes, parent::AbstractFrameAxes) where {T, O}

    # Get nutation - precession angles 
    p = Orient.PlanetsPrecessionNutation(Frames.point_id(center), data)

    # Build transformation
    # In this case the transformation is a static rotation matrix 
    dcm = _orient_bci2000(p)

    # Insert new axes in the frame system 
    add_axes_fixedoffset!(frames, axes, parent, dcm)
end