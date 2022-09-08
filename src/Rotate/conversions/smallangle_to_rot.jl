export smallangle_to_rot

"""
    smallangle_to_rot([T,] θx::Number, θy::Number, θz::Number[; normalize = true])

Create a rotation description of type `T` from three small rotations of angles
`θx`, `θy`, and `θz` [rad] about the axes X, Y, and Z, respectively.

The type `T` of the rotation description can be `DCM` or `Quaternion`. If the
type `T` is not specified, then if defaults to `DCM`.

If `T` is `DCM`, then the resulting matrix will be orthonormalized using the
`orthonormalize` function if the keyword `normalize` is `true`.
"""
@inline function smallangle_to_rot(
    θx::Number,
    θy::Number,
    θz::Number;
    normalize = true
)
    return smallangle_to_dcm(θx, θy, θz; normalize = normalize)
end

@inline function smallangle_to_rot(
    ::Type{DCM},
    θx::Number,
    θy::Number,
    θz::Number;
    normalize = true
)
    return smallangle_to_dcm(θx, θy, θz; normalize = normalize)
end

@inline function smallangle_to_rot(
    ::Type{Quaternion},
    θx::Number,
    θy::Number,
    θz::Number
)
    return smallangle_to_quat(θx, θy, θz)
end
