
"""
    add_axes_bcrtod!(fr, name, id, center; deriv=true)

Add Body-Centered Rotating (BCR), True-of-Date (TOD) axes with `name` and `id` to `fr`. 
The center point (i.e., the reference body) is `center`.

These axes are the equivalent of SPICE's `IAU_<BODY_NAME>` frames.

!!! warning 
    The parent axes are automatically set to the ICRF (ID = $(AXESID_ICRF)). If the 
    ICRF is not defined in `fr`, an error is thrown.

### See also 
See also [`add_axes_rotating!`](@ref), [`add_axes_bci2000!`](@ref).
"""
function add_axes_bcrtod!(
    fr::FrameSystem, name::Symbol, id::Int, center; deriv::Bool=false
)
    # create val dispatch
    cid = point_id(fr, center)
    vid = Val(cid)

    if length(methods(body_rotational_elements, [Number, Val{cid}])) != 1
        throw(
            MethodError(
                body_rotational_elements,
                "there must be a unique method defined for $cid"
            )
        )
    end

    if !(has_axes(fr, AXESID_ICRF))
        throw(
            ErrorException(
                "Body-Centered Inertial (BCI) at J2000 axes can only be defined" *
                " w.r.t. the ICRF (ID = $(AXESID_ICRF)), which is not defined in" *
                " the current frames graph."
            )
        )
    end

    if !deriv
        add_axes_rotating!(fr, name, id, AXESID_ICRF, t -> _bcrtod(t, vid))
    else
        add_axes_rotating!(
            fr, name, id, AXESID_ICRF,
            t -> _bcrtod(t, vid), t -> _∂bcrtod(t, vid), t -> _∂²bcrtod(t, vid), t -> _∂³bcrtod(t, vid)
        )
    end
end

"""
    add_axes_bci2000!(fr, axes::AbstractFrameAxes, center, data)

Add Body-Centered Inertial (BCI) axes at J2000 with `name` and `id` to `fr`. 
The center point (i.e., the reference body) is `center`.

!!! warning 
    The parent axes are automatically set to the ICRF (ID = $(AXESID_ICRF)). If the 
    ICRF is not defined in `fr`, an error is thrown.
    
### See also 
See also [`add_axes_fixedoffset!`](@ref), [`add_axes_bcrtod!`](@ref).
"""
function add_axes_bci2000!(fr::FrameSystem, name::Symbol, id::Int, center)

    # create val dispatch
    cid = point_id(fr, center)
    vid = Val(cid)

    if length(methods(body_rotational_elements, [Number, Val{cid}])) != 1
        throw(
            ErrorException(
                "there must be a unique method defined for $cid"
            )
        )
    end

    if !(has_axes(fr, AXESID_ICRF))
        throw(
            ErrorException(
                "Body-Centered Inertial (BCI) at J2000 axes can only be defined" *
                " w.r.t. the ICRF (ID = $(AXESID_ICRF)), which is not defined in" *
                " the current frames graph."
            )
        )
    end

    # evaluate rotational elements and build the DCM 
    α2000, δ2000, _ = body_rotational_elements(0, vid)
    dcm = angle_to_dcm(π / 2 + α2000, π / 2 - δ2000, :ZX)

    # insert the new axes
    return add_axes_fixedoffset!(fr, name, id, AXESID_ICRF, dcm)

end

function _bcrtod(seconds, val)
    α, δ, w = body_rotational_elements(seconds / CENTURY2SEC, val)
    return angle_to_dcm(π / 2 + α, π / 2 - δ, w, :ZXZ)
end

function _∂bcrtod(seconds, val)
    T = seconds / Tempo.CENTURY2SEC
    α, δ, w = body_rotational_elements(T, val)
    dα, dδ, dw = ∂body_rotational_elements(T, val)

    R = angle_to_dcm(π / 2 + α, π / 2 - δ, w, :ZXZ)
    dR = angle_to_δdcm((π / 2 + α, dα), (π / 2 - δ, -dδ), (w, dw), :ZXZ)
    return R, dR
end

function _∂²bcrtod(seconds, val)
    T = seconds / Tempo.CENTURY2SEC
    α, δ, w = body_rotational_elements(T, val)
    dα, dδ, dw = ∂body_rotational_elements(T, val)
    d²α, d²δ, d²w = ∂²body_rotational_elements(T, val)

    R = angle_to_dcm(π / 2 + α, π / 2 - δ, w, :ZXZ)
    dR = angle_to_δdcm((π / 2 + α, dα), (π / 2 - δ, -dδ), (w, dw), :ZXZ)
    d²R = angle_to_δ²dcm((π / 2 + α, dα, d²α), (π / 2 - δ, -dδ, -d²δ), (w, dw, d²w), :ZXZ)
    return R, dR, d²R
end

function _∂³bcrtod(seconds, val)
    T = seconds / Tempo.CENTURY2SEC
    α, δ, w = body_rotational_elements(T, val)
    dα, dδ, dw = ∂body_rotational_elements(T, val)
    d²α, d²δ, d²w = ∂²body_rotational_elements(T, val)
    d³α, d³δ, d³w = ∂³body_rotational_elements(T, val)

    R = angle_to_dcm(π / 2 + α, π / 2 - δ, w, :ZXZ)
    dR = angle_to_δdcm((π / 2 + α, dα), (π / 2 - δ, -dδ), (w, dw), :ZXZ)
    d²R = angle_to_δ²dcm((π / 2 + α, dα, d²α), (π / 2 - δ, -dδ, -d²δ), (w, dw, d²w), :ZXZ)
    d³R = angle_to_δ³dcm((π / 2 + α, dα, d²α, d³α), (π / 2 - δ, -dδ, -d²δ, -d³δ), (w, dw, d²w, d³w), :ZXZ)
    return R, dR, d²R, d³R
end