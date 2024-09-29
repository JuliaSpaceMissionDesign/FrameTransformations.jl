
# ------------------------------------------------------------------------------------------
# PROMOTIONS
# ------------------------------------------------------------------------------------------

# Returns the inner datatype of a given DCM 
dcm_eltype(::Union{DCM{T},Type{DCM{T}}}) where {T} = T

# Returns a promoted type for a given tuple of DCMs 
@generated function promote_dcm_eltype(::Union{T,Type{T}}) where {T<:Tuple}
    t = Union{}
    for i in 1:length(T.parameters)
        tmp = dcm_eltype(Base.unwrapva(T.parameters[i]))
        t = :(promote_type($t, $tmp))
    end
    return quote
        Base.@_inline_meta
        $t
    end
end

# ------------------------------------------------------------------------------------------
# TYPE DEF
# ------------------------------------------------------------------------------------------

struct Rotation{S,T} <: AbstractArray{T,1}
    m::NTuple{S,DCM{T}}

    function Rotation(tup::NTuple{S,Any}) where {S}
        T = promote_dcm_eltype(tup)
        return new{S,T}(tup)
    end
end

@inline order(::Rotation{S,<:Any}) where {S} = S

# Julia API
Base.size(::Rotation{S,<:Any}) where {S} = (S,)
Base.getindex(R::Rotation, i) = R.m[i]
Base.length(::Rotation{S}) where {S} = S

# ------------------------------------------------------------------------------------------
# CONSTRUCTORS
# ------------------------------------------------------------------------------------------

# Varargs constructor
function Rotation(args::Vararg{Any,S}) where {S}
    return Rotation(args)
end

# Constructor with filter and auto-fill of missing DCMS 
@generated function Rotation{S}(dcms::DCM...) where {S}
    D = length(dcms)
    T = Expr(:call, :promote_dcm_eltype, :dcms)
    expr = Expr(:call, :tuple)
    for i in 1:min(S, D)
        push!(expr.args, Expr(:ref, :dcms, i))
    end
    for _ in 1:(S-D)
        push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, :zero, T), :I)))
    end
    return quote
        @inbounds Rotation($(expr))
    end
end

@generated function Rotation{S1}(dcms::NTuple{S2,DCM{T}}) where {S1,S2,T}
    expr = Expr(:call, :tuple)
    for i in 1:min(S1, S2)
        push!(expr.args, Expr(:ref, :dcms, i))
    end
    for _ in 1:(S1-S2)
        push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, :zero, T), :I)))
    end
    return quote
        @inbounds Rotation($(expr))
    end
end

# Constructor for S-order identity rotations! 
@generated function Rotation{S}(::UniformScaling{T}) where {S,T}
    expr = Expr(:call, :tuple)
    for i in 1:S
        if i == 1
            push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, :one, T), :I)))
        else
            push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, :zero, T), :I)))
        end
    end
    return quote
        @inbounds Rotation($(expr))
    end
end

@generated function Rotation{S,T}(::UniformScaling) where {S,T}
    expr = Expr(:call, :tuple)
    for i in 1:S
        if i == 1
            push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, T, 1), :I)))
        else
            push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, T, 0), :I)))
        end
    end
    return quote
        @inbounds Rotation($(expr))
    end
end

# Convert a Rotation to a different order
@generated function Rotation{S1}(rot::Rotation{S2,T}) where {S1,S2,T}
    expr = Expr(:call, :tuple)
    for i in 1:min(S1, S2)
        push!(expr.args, Expr(:ref, :rot, i))
    end
    for _ in 1:(S1-S2)
        push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, :zero, T), :I)))
    end
    return quote
        @inbounds Rotation($(expr))
    end
end

# Convert a rotation to a different order and type
@generated function Rotation{S1,T}(rot::Rotation{S2}) where {S1,S2,T}

    expr = Expr(:call, :tuple)
    for i in 1:min(S1, S2)
        push!(expr.args, Expr(:., T, Expr(:tuple, Expr(:ref, :rot, i))))
    end
    for _ in 1:(S1-S2)
        push!(expr.args, Expr(:call, :DCM, Expr(:call, :(*), Expr(:call, :zero, T), :I)))
    end
    return quote
        @inbounds Rotation($(expr))
    end
end

function Rotation(m::DCM{T}, ω::AbstractVector) where {T}
    dm = DCM(ddcm(m, SVector(ω)))
    return Rotation((m, dm))
end

@inline Rotation{S}(rot::Rotation{S}) where {S} = rot

# ------------------------------------------------------------------------------------------
# TYPE CONVERSIONS
# ------------------------------------------------------------------------------------------

# Convert a Rotation to a tuple 
@generated function Base.Tuple(rot::Rotation{S,T}) where {S,T}

    expr = Expr(:call, :tuple)
    for j in 1:(3S)
        Oⱼ = (j - 1) ÷ 3 + 1
        for i in 1:(3S)
            Oᵢ = (i - 1) ÷ 3 + 1
            if Oⱼ > Oᵢ
                push!(expr.args, Expr(:call, :zero, T))
            else
                row = i - 3 * (Oᵢ - 1)
                col = j - 3 * (Oⱼ - 1)
                rom = Oᵢ - Oⱼ + 1
                push!(
                    expr.args,
                    Expr(
                        :ref,
                        Expr(:ref, :rot, rom), row, col
                    )
                )
            end
        end
    end

    return quote
        Base.@_inline_meta
        @inbounds $(expr)
    end
end

# Generic Rotation-to-StaticArrays conversions
@inline function (::Type{SA})(rot::Rotation) where {SA<:StaticArray}
    return SA(Tuple(rot))
end

# ------------------------------------------------------------------------------------------
# OPERATIONS
# ------------------------------------------------------------------------------------------

# Inverse
Base.inv(rot::Rotation) = _inverse_rotation(rot)

@generated function _inverse_rotation(rot::Rotation{S,T}) where {S,T}
    expr = Expr(:call, :Rotation,)
    for i in 1:S
        push!(
            expr.args, Expr(:call, :adjoint, Expr(:ref, :rot, i))
        )
    end
    return quote
        @inbounds $(expr)
    end
end

# Product between Rotations
@inline Base.:*(r1::Rotation{S,<:Any}, r2::Rotation{S,<:Any}) where {S} = _compose_rotation(r1, r2)
function Base.:*(::Rotation{S1,<:Any}, ::Rotation{S2,<:Any}) where {S1,S2}
    throw(DimensionMismatch("Cannot multiply two `Rotation` types of order $S1 and $S2"))
end

@generated function _compose_rotation(A::Rotation{S,<:Any}, B::Rotation{S,<:Any}) where {S}
    expr = Expr(:call, :Rotation)
    for i in 1:S
        sum_expr = Expr(:call, :+)
        for j in 1:i
            c = binomial(i - 1, j - 1)
            aᵢ = Expr(:ref, :A, i - j + 1)
            bᵢ = Expr(:ref, :B, j)
            push!(sum_expr.args, Expr(:call, :*, c, aᵢ, bᵢ))
        end
        push!(expr.args, sum_expr)
    end

    return quote
        @inbounds $(expr)
    end
end

@inline Base.:*(A::Rotation{S}, v::Translation{S}) where {S} = _apply_rotation(A, v)
@inline function Base.:*(::Rotation{S1}, ::Translation{S2}) where {S1,S2}
    throw(DimensionMismatch("Cannot apply Rotation of order $S1 to Translation of order $S2"))
end

# Product between Rotation and a Translation
@generated function _apply_rotation(R::Rotation{S,Nr}, v::Translation{S,Nv}) where {S,Nr,Nv}
    # Apply rotation on a translation vector with the same size
    #
    #                    n
    # z(n) = R(n)*v(1) + ∑ binomial(n, k) * R(n-k) * v(k)
    #                   k=1

    expr = Expr(:call, :tuple)
    for i in 1:S
        sumexpr = Expr(:call, :+)
        push!(expr.args, sumexpr)
        push!(sumexpr.args, Expr(:call, :*, Expr(:ref, :R, i), Expr(:ref, :v, 1)))

        for j in 1:i-1
            push!(
                sumexpr.args,
                Expr(
                    :call, :*,
                    binomial(i - 1, j),
                    Expr(:ref, :R, i - j),
                    Expr(:ref, :v, j + 1),
                )
            )
        end
    end
    return quote
        Base.@_inline_meta
        @inbounds Translation($(expr))
    end
end

@inline Base.:*(A::Rotation, v::SVector) = _apply_rotation(A, v)

# Compute product between Rotation and a "proper" SVector (returning a Translation)
@generated function _apply_rotation(R::Rotation{Sr,Nr}, v::SVector{Sv,Nv}) where {Sr,Sv,Nr,Nv}

    if Sv != 3Sr
        throw(
            DimensionMismatch(
                "Cannot apply Rotation of order $Sr to a $(Sv) vector",
            )
        )
    end

    expr = Expr(
        :call,
        Expr(:curly, :SVector, Sv, Nv),
        Expr(
            :call,
            :_apply_rotation, :R,
            Expr(:call, Expr(:curly, :Translation, Sr), :v)
        )
    )

    return quote
        Base.@_inline_meta
        @inbounds $expr
    end
end