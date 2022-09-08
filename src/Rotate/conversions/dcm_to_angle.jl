export dcm_to_angle

"""
    dcm_to_angle(dcm::DCM, seq::Symbol=:ZYX)

Convert the `dcm` to Euler Angles (see [`EulerAngles`](@ref)) given a rotation
sequence `seq`.

The rotation sequence is defined by a `:Symbol`. The possible values are:
`:XYX`, `XYZ`, `:XZX`, `:XZY`, `:YXY`, `:YXZ`, `:YZX`, `:YZY`, `:ZXY`, `:ZXZ`,
`:ZYX`, and `:ZYZ`. If no value is specified, then it defaults to `:ZYX`.

# Gimbal-lock and special cases

If the rotations are about three different axes, *e.g.* `:XYZ`, `:ZYX`, etc.,
then a second rotation of `±90˚` yields a gimbal-lock. This means that the
rotations between the first and third axes have the same effect. In this case,
the net rotation angle is assigned to the first rotation, and the angle of the
third rotation is set to 0.

If the rotations are about two different axes, *e.g.* `:XYX`, `:YXY`, etc., then
a rotation about the duplicated axis yields multiple representations. In this
case, the entire angle is assigned to the first rotation and the third rotation
is set to 0.
"""
function dcm_to_angle(dcm::DCM{T}, seq::Symbol=:ZYX) where T<:Number
    if seq == :ZYX
        # Check for singularities.
        if !(abs(dcm[1, 3]) ≥ 1 - eps())
            return EulerAngles(
                _mod_atan(+dcm[1, 2], +dcm[1, 1]),
                asin(-dcm[1, 3]),
                _mod_atan(+dcm[2, 3], +dcm[3, 3]),
                seq
            )
        else
            return EulerAngles(
                _mod_atan(-dcm[2, 1], +dcm[2, 2]),
                _mod_asin(-dcm[1, 3]),
                T(0),
                seq
            )
        end
    elseif seq == :XYX
        # Check for singularities.
        if !(abs(dcm[1, 1]) ≥ 1 - eps())
            return EulerAngles(
                _mod_atan(+dcm[1, 2], -dcm[1, 3]),
                acos(+dcm[1, 1]),
                _mod_atan(+dcm[2, 1], +dcm[3, 1]),
                seq
            )
        else
            return EulerAngles(
                _mod_atan(+dcm[2, 3], +dcm[2, 2]),
                _mod_acos(dcm[1, 1]),
                T(0),
                seq
            )
        end
    elseif seq == :XYZ
        # Check for singularities.
        if !(abs(dcm[3, 1]) ≥ 1 - eps())
            return EulerAngles(
                _mod_atan(-dcm[3, 2], +dcm[3, 3]),
                asin(+dcm[3, 1]),
                _mod_atan(-dcm[2, 1], +dcm[1, 1]),
                seq
            )
        else
            return EulerAngles(
                _mod_atan(+dcm[2, 3], +dcm[2, 2]),
                _mod_asin(+dcm[3, 1]),
                T(0),
                seq
            )
        end
    elseif seq == :XZX
        # Check for singularities.
        if !(abs(dcm[1, 1]) ≥ 1 - eps())
            return EulerAngles(
                _mod_atan(+dcm[1, 3], +dcm[1, 2]),
                acos(+dcm[1, 1]),
                _mod_atan(+dcm[3, 1], -dcm[2, 1]),
                seq
            )
        else
            return EulerAngles(
                _mod_atan(-dcm[3, 2], +dcm[3, 3]),
                _mod_acos(dcm[1, 1]),
                T(0),
                seq
            )
        end
    elseif seq == :XZY
        # Check for singularities.
        if !(abs(dcm[2, 1]) ≥ 1 - eps())
            return EulerAngles(
                _mod_atan(+dcm[2, 3], +dcm[2, 2]),
                asin(-dcm[2, 1]),
                _mod_atan(+dcm[3, 1], +dcm[1, 1]),
                seq
            )
        else
            return EulerAngles(
                _mod_atan(-dcm[3, 2], +dcm[3, 3]),
                _mod_asin(-dcm[2, 1]),
                T(0),
                seq
            )
        end
    elseif seq == :YXY
        # Check for singularities.
        if !(abs(dcm[2, 2]) ≥ 1 - eps())
            return EulerAngles(
                _mod_atan(+dcm[2, 1], +dcm[2, 3]),
                acos(+dcm[2, 2]),
                _mod_atan(+dcm[1, 2], -dcm[3, 2]),
                seq
            )
        else
            return EulerAngles(
                _mod_atan(-dcm[1, 3], +dcm[1, 1]),
                _mod_acos(dcm[2, 2]),
                T(0),
                seq
            )
        end
    elseif seq == :YXZ
        if !(abs(dcm[3, 2]) ≥ 1 - eps())
            return EulerAngles(
                _mod_atan(+dcm[3, 1], +dcm[3, 3]),
                asin(-dcm[3, 2]),
                _mod_atan(+dcm[1, 2], +dcm[2, 2]),
                seq
            )
        else
            return EulerAngles(
                _mod_atan(-dcm[1, 3], +dcm[1, 1]),
                _mod_asin(-dcm[3, 2]),
                T(0),
                seq
            )
        end
    elseif seq == :YZX
        # Check for singularities.
        if !(abs(dcm[1, 2]) ≥ 1 - eps())
            return EulerAngles(
                _mod_atan(-dcm[1, 3], +dcm[1, 1]),
                asin(+dcm[1, 2]),
                _mod_atan(-dcm[3, 2], +dcm[2, 2]),
                seq
            )
        else
            return EulerAngles(
                _mod_atan(+dcm[3, 1], +dcm[3, 3]),
                _mod_asin(+dcm[1, 2]),
                T(0),
                seq
            )
        end
    elseif seq == :YZY
        # Check for singularities.
        if !(abs(dcm[2, 2]) ≥ 1 - eps())
            return EulerAngles(
                _mod_atan(+dcm[2, 3], -dcm[2, 1]),
                acos(+dcm[2, 2]),
                _mod_atan(+dcm[3, 2], +dcm[1, 2]),
                seq
            )
        else
            return EulerAngles(
                _mod_atan(dcm[3, 1], dcm[3, 3]),
                _mod_acos(dcm[2, 2]),
                T(0),
                seq
            )
        end
    elseif seq == :ZXY
        # Check for singularities.
        if !(abs(dcm[2, 3]) ≥ 1 - eps())
            return EulerAngles(
                _mod_atan(-dcm[2, 1], +dcm[2, 2]),
                asin(+dcm[2, 3]),
                _mod_atan(-dcm[1, 3], +dcm[3, 3]),
                seq
            )
        else
            return EulerAngles(
                _mod_atan(+dcm[1, 2], +dcm[1, 1]),
                _mod_asin(+dcm[2, 3]),
                T(0),
                seq
            )
        end
    elseif seq == :ZXZ
        # Check for singularities.
        if !(abs(dcm[3, 3]) ≥ 1 - eps())
            return EulerAngles(
                _mod_atan(+dcm[3, 1], -dcm[3, 2]),
                acos(+dcm[3, 3]),
                _mod_atan(+dcm[1, 3], +dcm[2, 3]),
                seq
            )
        else
            return EulerAngles(
                _mod_atan(dcm[1, 2], dcm[1, 1]),
                _mod_acos(dcm[3, 3]),
                T(0),
                seq
            )
        end
    elseif seq == :ZYZ
        # Check for singularities.
        if !(abs(dcm[3, 3]) ≥ 1 - eps())
            return EulerAngles(
                _mod_atan(+dcm[3, 2], +dcm[3, 1]),
                acos(+dcm[3, 3]),
                _mod_atan(+dcm[2, 3], -dcm[1, 3]),
                seq
            )
        else
            return EulerAngles(
                _mod_atan(-dcm[2, 1], dcm[2, 2]),
                _mod_acos(dcm[3, 3]),
                T(0),
                seq
            )
        end
    else
        throw(ArgumentError("The rotation sequence :$seq is not valid."))
    end
end

################################################################################
#                              Private functions
################################################################################

# This modified function computes exactly what `atan(y,x)` computes except that
# it will neglect signed zeros. Hence:
#
#   _mod_atan(0.0, -0.0) = _mod_atan(-0.0, 0.0) = 0.0
#
# The signed zero can lead to problems when converting from DCM to Euler angles.
_mod_atan(y::T, x::T) where T<:Number = atan(y + T(0), x + T(0))

# This modified function computes the `acos(x)` if `|x| <= 1` and computes
# `acos( sign(x) )`  if `|x| > 1` to avoid numerical errors when converting DCM
# to  Euler Angles.
function _mod_acos(x::T) where T<:Number
    if x > 1
        return float(T(0))
    elseif x < -1
        return float(T(π))
    else
        return acos(x)
    end
end

# This modified function computes the `asin(x)` if `|x| <= 1` and computes
# `asin( sign(x) )`  if `|x| > 1` to avoid numerical errors when converting DCM
# to  Euler Angles.
function _mod_asin(x::T) where T<:Number
    if x > 1
        return +float(T(π / 2))
    elseif x < -1
        return -float(T(π / 2))
    else
        return asin(x)
    end
end
