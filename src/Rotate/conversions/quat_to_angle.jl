export quat_to_angle

"""
    quat_to_angle(q::Quaternion, seq::Symbol = :ZYX)

Convert the quaternion `q` to Euler Angles (see [`EulerAngles`](@ref)) given a
rotation sequence `seq`.

The rotation sequence is defined by a `:Symbol`. The possible values are:
`:XYX`, `XYZ`, `:XZX`, `:XZY`, `:YXY`, `:YXZ`, `:YZX`, `:YZY`, `:ZXY`, `:ZXZ`,
`:ZYX`, and `:ZYZ`. If no value is specified, then it defaults to `:ZYX`.
"""
function quat_to_angle(q::Quaternion, seq::Symbol = :ZYX)
    # Convert the quaternion to DCM.
    dcm = quat_to_dcm(q)

    # Convert the DCM to the Euler Angles.
    return dcm_to_angle(dcm, seq)
end
