export angle_to_quat

"""
    angle_to_quat(θ₁::T1[, θ₂::T2[, θ₃::T3]], seq::Symbol = :ZYX) where {T1<:Number, T2<:Number, T3<:Number}
    angle_to_quat(eulerang::EulerAngles)

Create a quaternion that perform a set of rotations (`θ₁`, `θ₂`, `θ₃`) about the
coordinate axes specified in `seq`.

The input values of the origin Euler angles can also be passed inside the
structure `Θ` (see [`EulerAngles`](@ref)).

The rotation sequence is defined by a `Symbol` specifing the rotation axes. The
possible values depends on the number of rotations as follows:

- **1 rotation** (`θ₁`): `:X`, `:Y`, or `:Z`.
- **2 rotations** (`θ₁`, `θ₂`): `:XY`, `:XZ`, `:YX`, `:YZ`, `:ZX`, or `:ZY`.
- **3 rotations** (`θ₁`, `θ₂`, `θ₃`): `:XYX`, `XYZ`, `:XZX`, `:XZY`, `:YXY`,
    `:YXZ`, `:YZX`, `:YZY`, `:ZXY`, `:ZXZ`, `:ZYX`, or `:ZYZ`

!!! note
    The type of the new quaternion will be obtained by promiting `T1`, `T2`, and
    `T3`.

# Remarks

This function assigns `q = q1 * q2 * q3` in which `qi` is the quaternion related
with the *i*-th rotation, `i Є [1,2,3]`. If the *i*-th rotation is not specified,
then `qi = Quaternion(I)`.
"""
function angle_to_quat(θ::Number, seq::Symbol)
    # Compute the sines and cosines of half angle.
    s, c = sincos(θ / 2)

    # Make sure that the real part is always positive.
    if c < 0
        c = -c
        s = -s
    end

    if seq == :X
        return Quaternion(c, s, 0, 0)
    elseif seq == :Y
        return Quaternion(c, 0, s, 0)
    elseif seq == :Z
        return Quaternion(c, 0, 0, s)
    else
        throw(ArgumentError("seq must be :X, :Y, or :Z"))
    end
end

function angle_to_quat(
    θ₁::T1,
    θ₂::T2,
    seq::Symbol
) where {T1<:Number, T2<:Number}
    T = promote_type(T1, T2)

    # Compute the sines and cosines of half angle.
    s₁, c₁ = sincos(T(θ₁) / 2)
    s₂, c₂ = sincos(T(θ₂) / 2)

    # When we have two rotations, the `q0` component is always the same.
    q0 = c₁ * c₂

    s = (q0 < 0) ? -1 : +1

    if seq == :XY
        return Quaternion(
            s * q0,
            s * (s₁ * c₂),
            s * (c₁ * s₂),
            s * (s₁ * s₂)
        )
    elseif seq == :XZ
        return Quaternion(
            s * q0,
            s * ( s₁ * c₂),
            s * (-s₁ * s₂),
            s * ( c₁ * s₂)
        )
    elseif seq == :YX
        return Quaternion(
            s * q0,
            s * ( c₁ * s₂),
            s * ( s₁ * c₂),
            s * (-s₁ * s₂)
        )
    elseif seq == :YZ
        return Quaternion(
            s * q0,
            s * (s₁ * s₂),
            s * (s₁ * c₂),
            s * (c₁ * s₂)
        )
    elseif seq == :ZX
        return Quaternion(
            s * q0,
            s * (c₁ * s₂),
            s * (s₁ * s₂),
            s * (s₁ * c₂)
        )
    elseif seq == :ZY
        return Quaternion(
            s * q0,
            s * (-s₁ * s₂),
            s * ( c₁ * s₂),
            s * ( s₁ * c₂)
        )
    else
        throw(ArgumentError("The rotation sequence :$seq is not valid."))
    end
end

function angle_to_quat(
    θ₁::T1,
    θ₂::T2,
    θ₃::T3,
    seq::Symbol = :ZYX
) where {T1<:Number, T2<:Number, T3<:Number}
    T = promote_type(T1, T2, T3)

    # Compute the sines and cosines of half angle.
    s₁, c₁ = sincos(T(θ₁) / 2)
    s₂, c₂ = sincos(T(θ₂) / 2)
    s₃, c₃ = sincos(T(θ₃) / 2)

    if seq == :ZYX
        q0 = c₁ * c₂ * c₃ + s₁ * s₂ * s₃

        s = (q0 < 0) ? -1 : +1

        return Quaternion(
            s * q0,
            s * (c₁ * c₂ * s₃ - s₁ * s₂ * c₃),
            s * (c₁ * s₂ * c₃ + s₁ * c₂ * s₃),
            s * (s₁ * c₂ * c₃ - c₁ * s₂ * s₃)
        )
    elseif seq == :XYX
        q0 = c₁ * c₂ * c₃ - s₁ * c₂ * s₃

        s = (q0 < 0) ? -1 : +1

        return Quaternion(
            s * q0,
            s * (c₁ * c₂ * s₃ + s₁ * c₂ * c₃),
            s * (c₁ * s₂ * c₃ + s₁ * s₂ * s₃),
            s * (s₁ * s₂ * c₃ - c₁ * s₂ * s₃)
        )
    elseif seq == :XYZ
        q0 = c₁ * c₂ * c₃ - s₁ * s₂ * s₃

        s = (q0 < 0) ? -1 : +1

        return Quaternion(
            s * q0,
            s * (s₁ * c₂ * c₃ + c₁ * s₂ * s₃),
            s * (c₁ * s₂ * c₃ - s₁ * c₂ * s₃),
            s * (c₁ * c₂ * s₃ + s₁ * s₂ * c₃)
        )
    elseif seq == :XZX
        q0 = c₁ * c₂ * c₃ - s₁ * c₂ * s₃

        s = (q0 < 0) ? -1 : +1

        return Quaternion(
            s * q0,
            s * (c₁ * c₂ * s₃ + s₁ * c₂ * c₃),
            s * (c₁ * s₂ * s₃ - s₁ * s₂ * c₃),
            s * (c₁ * s₂ * c₃ + s₁ * s₂ * s₃)
        )
    elseif seq == :XZY
        q0 = c₁ * c₂ * c₃ + s₁ * s₂ * s₃

        s = (q0 < 0) ? -1 : +1

        return Quaternion(
            s * q0,
            s * (s₁ * c₂ * c₃ - c₁ * s₂ * s₃),
            s * (c₁ * c₂ * s₃ - s₁ * s₂ * c₃),
            s * (c₁ * s₂ * c₃ + s₁ * c₂ * s₃)
        )
    elseif seq == :YXY
        q0 = c₁ * c₂ * c₃ - s₁ * c₂ * s₃

        s = (q0 < 0) ? -1 : +1

        return Quaternion(
            s * q0,
            s * (c₁ * s₂ * c₃ + s₁ * s₂ * s₃),
            s * (c₁ * c₂ * s₃ + s₁ * c₂ * c₃),
            s * (c₁ * s₂ * s₃ - s₁ * s₂ * c₃)
        )
    elseif seq == :YXZ
        q0 = c₁ * c₂ * c₃ + s₁ * s₂ * s₃

        s = (q0 < 0) ? -1 : +1

        return Quaternion(
            s * q0,
            s * (c₁ * s₂ * c₃ + s₁ * c₂ * s₃),
            s * (s₁ * c₂ * c₃ - c₁ * s₂ * s₃),
            s * (c₁ * c₂ * s₃ - s₁ * s₂ * c₃)
        )
    elseif seq == :YZX
        q0 = c₁ * c₂ * c₃ - s₁ * s₂ * s₃

        s = (q0 < 0) ? -1 : +1

        return Quaternion(
            s * q0,
            s * (c₁ * c₂ * s₃ + s₁ * s₂ * c₃),
            s * (s₁ * c₂ * c₃ + c₁ * s₂ * s₃),
            s * (c₁ * s₂ * c₃ - s₁ * c₂ * s₃)
        )
    elseif seq == :YZY
        q0 = c₁ * c₂ * c₃ - s₁ * c₂ * s₃

        s = (q0 < 0) ? -1 : +1

        return Quaternion(
            s * q0,
            s * (s₁ * s₂ * c₃ - c₁ * s₂ * s₃),
            s * (c₁ * c₂ * s₃ + s₁ * c₂ * c₃),
            s * (c₁ * s₂ * c₃ + s₁ * s₂ * s₃)
        )
    elseif seq == :ZXY
        q0 = c₁ * c₂ * c₃ - s₁ * s₂ * s₃

        s = (q0 < 0) ? -1 : +1

        return Quaternion(
            s * q0,
            s * (c₁ * s₂ * c₃ - s₁ * c₂ * s₃),
            s * (c₁ * c₂ * s₃ + s₁ * s₂ * c₃),
            s * (s₁ * c₂ * c₃ + c₁ * s₂ * s₃)
        )
    elseif seq == :ZXZ
        q0 = c₁ * c₂ * c₃ - s₁ * c₂ * s₃

        s = (q0 < 0) ? -1 : +1

        return Quaternion(
            s * q0,
            s * (c₁ * s₂ * c₃ + s₁ * s₂ * s₃),
            s * (s₁ * s₂ * c₃ - c₁ * s₂ * s₃),
            s * (c₁ * c₂ * s₃ + s₁ * c₂ * c₃)
        )
    elseif seq == :ZYZ
        q0 = c₁ * c₂ * c₃ - s₁ * c₂ * s₃

        s = (q0 < 0) ? -1 : +1

        return Quaternion(
            s * q0,
            s * (c₁ * s₂ * s₃ - s₁ * s₂ * c₃),
            s * (c₁ * s₂ * c₃ + s₁ * s₂ * s₃),
            s * (c₁ * c₂ * s₃ + s₁ * c₂ * c₃)
        )
    else
        throw(ArgumentError("The rotation sequence :$seq is not valid."))
    end
end

angle_to_quat(Θ::EulerAngles) = angle_to_quat(Θ.a1, Θ.a2, Θ.a3, Θ.seq)
