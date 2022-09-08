export orthonormalize, δdcm

# StaticArrays API

@inline Base.@propagate_inbounds function getindex(dcm::DCM, i::Int)
    return dcm.data[i]
end

function Base.Tuple(dcm::DCM)
    return dcm.data
end

function similar_type(::Type{A}, ::Type{T}, s::Size{(3,3)}) where {A<:DCM, T}
    return DCM{T}
end

# Julia API

function summary(io::IO, ::DCM{T}) where T
    print(io, "DCM{" * string(T) * "}")
end

"""
    orthonormalize(dcm::DCM)

Perform the Gram-Schmidt orthonormalization process in the `dcm` and return the
new matrix.

!!! warning
    This function does not check if the columns of the input matrix span a
    three-dimensional space. If not, then the returned matrix should have `NaN`.
    Notice, however, that such input matrix is not a valid direction cosine
    matrix.
"""
function orthonormalize(dcm::DCM)
    e₁ = dcm[:, 1]
    e₂ = dcm[:, 2]
    e₃ = dcm[:, 3]

    en₁  = e₁ / norm(e₁)
    enj₂ = e₂ - (en₁ ⋅ e₂) * en₁
    en₂  = enj₂ / norm(enj₂)
    enj₃ = e₃ - (en₁ ⋅ e₃) * en₁
    enj₃ = enj₃ - (en₂ ⋅ enj₃) * en₂
    en₃  = enj₃ / norm(enj₃)

    return DCM(hcat(en₁, en₂, en₃))
end

"""
    δdcm(Dba::DCM, ω::AbstractArray)

Compute the time-derivative of the `dcm` that rotates a reference frame `a` into
alignment with the reference frame `b` in which the angular velocity of `b` with
respect to `a`, and represented in `b`, is `ω`.
"""
function δdcm(Dba::DCM, ω::AbstractArray)
    if length(ω) != 3 # Check the dimensions.
        throw(ArgumentError("The angular velocity vector must have three components."))
    end
    ωx = SMatrix{3, 3}(
          0  , -ω[3], +ω[2],
        +ω[3],   0  , -ω[1],
        -ω[2], +ω[1],   0,
    )'
    return -ωx * Dba
end
