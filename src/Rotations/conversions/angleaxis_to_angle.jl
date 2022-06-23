export angleaxis_to_angle

"""
    angleaxis_to_angle(θ::Number, v::AbstractVector, seq::Symbol)
    angleaxis_to_angle(av::EulerAngleAxis, seq::Symbol)

Convert the Euler angle `θ` [rad]  and Euler axis `v` to Euler angles with
rotation sequence `seq`.

Those values can also be passed inside the structure `av` (see
[`EulerAngleAxis`](@ref)).

The rotation sequence is defined by a `:Symbol`. The possible values are:
`:XYX`, `XYZ`, `:XZX`, `:XZY`, `:YXY`, `:YXZ`, `:YZX`, `:YZY`, `:ZXY`, `:ZXZ`,
`:ZYX`, and `:ZYZ`. If no value is specified, then it defaults to `:ZYX`.

!!! warning
    It is expected that the vector `v` is unitary. However, no verification is
    performed inside the function. The user must handle this situation.
"""
@inline function angleaxis_to_angle(θ::Number, v::AbstractVector, seq::Symbol)
    # Check the arguments.
    (length(v) ≠ 3) && throw(ArgumentError("The provided vector for the Euler axis must have 3 elements."))

    # First we convert to DCM then to Euler angles.
    return dcm_to_angle(angleaxis_to_dcm(θ, v), seq)
end

@inline function angleaxis_to_angle(av::EulerAngleAxis, seq::Symbol)
    return angleaxis_to_angle(av.a, av.v, seq)
end
