export angle_to_angleaxis

"""
    angle_to_angleaxis(θ₁::Number, θ₂::Number, θ₃::Number, seq::Symbol = :ZYX)
    angle_to_angleaxis(Θ::EulerAngles)

Convert the Euler angles `θ₁`, `θ₂`, and `θ₃` [rad] with the rotation sequence
`seq` to an Euler angle and axis representation.

Those values can also be passed inside the structure `Θ` (see
[`EulerAngles`](@ref)).

The rotation sequence is defined by a `:Symbol`. The possible values are:
`:XYX`, `XYZ`, `:XZX`, `:XZY`, `:YXY`, `:YXZ`, `:YZX`, `:YZY`, `:ZXY`, `:ZXZ`,
`:ZYX`, and `:ZYZ`. If no value is specified, then it defaults to `:ZYX`.

"""
@inline function angle_to_angleaxis(
    θ₁::Number,
    θ₂::Number,
    θ₃::Number,
    seq::Symbol = :ZYX
)
    # First convert to DCM and then to Euler angle and axis.
    return dcm_to_angleaxis(angle_to_dcm(θ₁, θ₂, θ₃, seq))
end

@inline function angle_to_angleaxis(Θ::EulerAngles)
    return angle_to_angleaxis(Θ.a1, Θ.a2, Θ.a3, Θ.seq)
end
