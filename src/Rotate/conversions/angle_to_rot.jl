export angle_to_rot

"""
    angle_to_rot([T,] θ₁::Number[, θ₂::Number[, θ₃::Number]], seq::Symbol)
    angle_to_rot([T,] Θ::EulerAngles)

Create a rotation description of type `T` that perform a set of rotations (`θ₁`,
`θ₂`, `θ₃`) about the coordinate axes specified in `seq`.

The input values of the origin Euler angles can also be passed inside the
structure `Θ` (see [`EulerAngles`](@ref)).

The rotation sequence is defined by a `Symbol` specifing the rotation axes. The
possible values depends on the number of rotations as follows:

- **1 rotation** (`θ₁`): `:X`, `:Y`, or `:Z`.
- **2 rotations** (`θ₁`, `θ₂`): `:XY`, `:XZ`, `:YX`, `:YZ`, `:ZX`, or `:ZY`.
- **3 rotations** (`θ₁`, `θ₂`, `θ₃`): `:XYX`, `XYZ`, `:XZX`, `:XZY`, `:YXY`,
    `:YXZ`, `:YZX`, `:YZY`, `:ZXY`, `:ZXZ`, `:ZYX`, or `:ZYZ`
"""
@inline function angle_to_rot(θ::Number, seq::Symbol)
    return angle_to_dcm(θ, seq)
end

@inline function angle_to_rot(::Type{DCM}, θ::Number, seq::Symbol)
    return angle_to_dcm(θ, seq)
end

@inline function angle_to_rot(::Type{Quaternion}, θ::Number, seq::Symbol)
    return angle_to_quat(θ, seq)
end

@inline function angle_to_rot(
    θ₁::Number,
    θ₂::Number,
    seq::Symbol
)
    return angle_to_dcm(θ₁, θ₂, seq)
end

@inline function angle_to_rot(
    ::Type{DCM},
    θ₁::Number,
    θ₂::Number,
    seq::Symbol
)
    return angle_to_dcm(θ₁, θ₂, seq)
end

@inline function angle_to_rot(
    ::Type{Quaternion},
    θ₁::Number,
    θ₂::Number,
    seq::Symbol
)
    return angle_to_quat(θ₁, θ₂, seq)
end

@inline function angle_to_rot(
    θ₁::Number,
    θ₂::Number,
    θ₃::Number,
    seq::Symbol
)
    return angle_to_dcm(θ₁, θ₂, θ₃, seq)
end

@inline function angle_to_rot(
    ::Type{DCM},
    θ₁::Number,
    θ₂::Number,
    θ₃::Number,
    seq::Symbol
)
    return angle_to_dcm(θ₁, θ₂, θ₃, seq)
end

@inline function angle_to_rot(
    ::Type{Quaternion},
    θ₁::Number,
    θ₂::Number,
    θ₃::Number,
    seq::Symbol
)
    return angle_to_quat(θ₁, θ₂, θ₃, seq)
end

@inline function angle_to_rot(Θ::EulerAngles)
    return angle_to_rot(Θ.a1, Θ.a2, Θ.a3, Θ.seq)
end

@inline function angle_to_rot(
    T::Union{Type{DCM}, Type{Quaternion}},
    Θ::EulerAngles
)
    return angle_to_rot(T, Θ.a1, Θ.a2, Θ.a3, Θ.seq)
end
