export inv_rotation

"""
    inv_rotation(R)

Compute the inverse rotation of `R`, which can be:

- A direction cosine matrix (`DCM`);
- An Euler angle and axis (`EulerAngleAxis`);
- A set of Euler anlges (`EulerAngles`); or
- A quaternion (`Quaternion`).

The output will have the same type as `R`.

!!! note
    If `R` is a DCM, than its transpose is computed instead of its inverse to
    reduce the computational burden. The both are equal if the DCM has unit
    norm. This must be verified by the user.

!!! note
    If `R` is a quaternion, than its conjugate is computed instead of its
    inverse to reduce the computational burden. The both are equal if the
    quaternion has unit norm. This must be verified by the used.
"""
@inline inv_rotation(D::DCM) = D'
@inline inv_rotation(ea::EulerAngleAxis) = inv(ea)
@inline inv_rotation(Θ::EulerAngles) = inv(Θ)
@inline inv_rotation(q::Quaternion) = conj(q)
