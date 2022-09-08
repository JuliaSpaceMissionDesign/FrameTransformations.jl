export angle_to_angle

"""
    angle_to_angle(θ₁::Number, θ₂::Number, θ₃::Number, seq_orig::Symbol, seq_dest::Symbol)
    angle_to_angle(Θ::EulerAngles, seq_dest::Symbol)

Convert the Euler angles `θ₁`, `θ₂`, and `θ₃` [rad] with the rotation sequence
`seq_orig` to a new set of Euler angles with rotation sequence
`seq_dest`.

The input values of the origin Euler angles can also be passed inside the
structure `Θ` (see [`EulerAngles`](@ref)).

The rotation sequence is defined by a `:Symbol`. The possible values are:
`:XYX`, `XYZ`, `:XZX`, `:XZY`, `:YXY`, `:YXZ`, `:YZX`, `:YZY`, `:ZXY`, `:ZXZ`,
`:ZYX`, and `:ZYZ`.

"""
@inline function angle_to_angle(
    θ₁::Number,
    θ₂::Number,
    θ₃::Number,
    seq_orig::Symbol,
    seq_dest::Symbol
)
    return dcm_to_angle(angle_to_dcm(θ₁, θ₂, θ₃, seq_orig), seq_dest)
end

@inline function angle_to_angle(Θ::EulerAngles, seq_dest::Symbol)
    return angle_to_angle(Θ.a1, Θ.a2, Θ.a3, Θ.seq, seq_dest)
end
