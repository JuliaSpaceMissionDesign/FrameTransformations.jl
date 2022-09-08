export compose_rotation

"""
    compose_rotation(R1, [, R2, R3, R4, R5, ...])

Compute a composed rotation using the rotations `R1`, `R2`, `R3`, `R4`, ..., in
the following order:

     First rotation
     |
     |
    R1 => R2 => R3 => R4 => ...
           |
           |
           Second rotation

The rotations can be described by:

- A direction cosine matrix ([`DCM`](@ref));
- An Euler angle and axis ([`EulerAngleAxis`](@ref));
- A set of Euler angles ([`EulerAngles`](@ref)); or
- A quaternion ([`Quaternion`](@ref)).

Notice, however, that all rotations **must be** of the same type (DCM or
quaternion).

The output will have the same type as the inputs.
"""
@inline compose_rotation(D::DCM) = D
@inline compose_rotation(D::DCM, Ds::DCM...) = compose_rotation(Ds...) * D

@inline compose_rotation(ea::EulerAngleAxis) = ea
@inline function compose_rotation(ea::EulerAngleAxis, eas::EulerAngleAxis...)
    return compose_rotation(eas...) * ea
end

@inline compose_rotation(Θ::EulerAngles) = Θ
@inline compose_rotation(Θ::EulerAngles, Θs::EulerAngles...) = compose_rotation(Θs...) * Θ

@inline compose_rotation(q::Quaternion) = q
@inline compose_rotation(q::Quaternion, qs::Quaternion...) = q * compose_rotation(qs...)

# This algorithm was proposed by @Per in
#
#   https://discourse.julialang.org/t/improve-the-performance-of-multiplication-of-an-arbitrary-number-of-matrices/10835/24

# Operator: ∘
# ==============================================================================

function ∘(R2::T1, R1::T2) where {
    T1<:Rotations,
    T2<:Rotations
}
    R1c = convert(T1, R1)
    return compose_rotation(R1c, R2)
end
