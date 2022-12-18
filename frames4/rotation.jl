using ReferenceFrameRotations
using StaticArrays
using LinearAlgebra: matprod

# -------------------------------------
# TYPES
# -------------------------------------

struct Rotation{S, N}
    m::NTuple{S, DCM{N}}
end

@generated function Rotation(dcms::DCM{N}...) where N
    S = length(dcms)
    expr = :(tuple($([Expr(:ref, :dcms, i) for i in 1:S]...)))
    return quote 
        @inbounds Rotation{$S, $N}($(expr))
    end
end

Base.getindex(R::Rotation, i) = R.m[i]

function Rotation(m::DCM{N}, ω::AbstractVector) where N
    dm = ddcm(m, SVector(ω))
    return Rotation{2, N}((m, dm))
end

# -------------------------------------
# OPERATIONS 
# -------------------------------------

# ---
# Inverse rotation
Base.inv(rot::Rotation) = _inverse_rot(rot)
@generated function _inverse_rot(rot::Rotation{S, N}) where {S, N}
    quote 
        @inbounds Rotation($((:(rot.m[$i]') for i in 1:S)...))
    end
end

# ---
# Compose rotations 
@inline Base.:*(A::Rotation{S, N}, B::Rotation{S, N}) where {S, N} = _compose_rot(A, B)

function Base.:*(A::Rotation{S1, N}, B::Rotation{S2, N}) where {S1, S2, N}
    throw(
        ArgumentError("Cannot multiply two `Rotation` of different order!")
    )
end

@generated function _compose_rot(A::Rotation{S, N}, B::Rotation{S, N}) where {S, N}
    expr = Expr(:call, :Rotation)

    for i in 1:S
        sum_expr = Expr(:call, :+, )
        for j in 1:i 
            c = binomial(i-1, j-1)
            ai = Expr(:ref, Expr(:., :A, QuoteNode(:m)), i-j+1)
            bi = Expr(:ref, Expr(:., :B, QuoteNode(:m)), j)

            push!(sum_expr.args, Expr(:call, :*, c, ai, bi))
        end
        push!(expr.args, sum_expr)    
    end

    return quote 
        @inbounds $(expr)
    end
end

# ---
# Apply rotation

@inline Base.:*(A::Rotation, b::AbstractVector) = _apply_rot(A, b)

# Function to compute product between Rotation and a generic vector
@generated function _apply_rot(A::Rotation{S, Na}, b::AbstractVector{Nb}) where {S, Na, Nb}

    exprs = [[Expr(:call, :+) for _ = 1:3] for _ = 1:S]

    for i = 1:S 
        for j = 1:i
            for k = 1:3 
                mi = i-j+1
                push!(exprs[i][k].args, 
                    StaticArrays.combine_products([:(A[$mi][$k, $w]*b[3*($j-1)+$w]) for w = 1:3]))
            end 
        end
    end

    sa = 3*S 
    retexpr = :(@inbounds return similar_type(b, T, Size($sa))(tuple($((exprs...)...))))
    return quote 
        Base.@_inline_meta
        length(b) != $sa && throw(DimensionMismatch(
            "Cannot multiply `Rotation` of size ($($sa), $($sa)) and a $(size(b)) vector"))

        T = Base.promote_op(matprod, Na, Nb)
        $retexpr
    end
end

# Function to compute product between Rotation and a static vector
@generated function _apply_rot(A::Rotation{S, Na}, b::StaticVector{sb, Nb}) where {S, sb, Na, Nb}
    sa = 3*S 
    sa != sb && throw(
        DimensionMismatch("Cannot multiply `Rotation` of size ($sa, $sa) and a $sb vector"))

    exprs = [[Expr(:call, :+) for _ = 1:3] for _ = 1:S]
    for i = 1:S 
        for j = 1:i
            for k = 1:3 
                push!(exprs[i][k].args, 
                    StaticArrays.combine_products([:(A[$i-$j+1][$k, $w]*b[3*($j-1)+$w]) for w = 1:3]))
            end 
        end
    end

    retexpr = :(@inbounds return similar_type(b, T, Size($sa))(tuple($((exprs...)...))))
    return quote 
        Base.@_inline_meta
        T = Base.promote_op(matprod, Na, Nb)
        $retexpr
    end
end
