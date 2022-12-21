using ReferenceFrameRotations
using StaticArrays
using LinearAlgebra: matprod, UniformScaling

import StaticArrays: similar_type, Size, MMatrix, SMatrix

# -------------------------------------
# TYPES
# -------------------------------------

struct Rotation{S, N}
    m::NTuple{S, DCM{N}}
    
    function Rotation(x::NTuple{S, Any}) where S 
        T = promote_dcm_eltype(x)
        return new{S, T}(x)
    end
end

order(::Rotation{S, <:Any}) where S = S

# Julia APIs 
Base.size(::Rotation{S, <:Any}) where S = (3S, 3S)
Base.getindex(R::Rotation, i) = R.m[i]

function Base.summary(io::IO, ::Rotation{S, N}) where {S, N}
    print(io, "Rotation{$S, $N}")
end

# ---
# Constructors 

# Generic rotation constructor 
@generated function Rotation(dcms::DCM...)
    S = length(dcms)
    expr = :(tuple($([Expr(:ref, :dcms, i) for i in 1:S]...)))
    return quote 
        @inbounds Rotation($(expr))
    end
end

# Constructor for S-order identity rotations! 
@generated function Rotation{S}(::UniformScaling{N}) where {S, N}
    expr = :(tuple($([i == 1 ? DCM(N(1)I) : DCM(N(0)I) for i in 1:S]...)))
    return quote 
        @inbounds Rotation($(expr))
    end
end

# Convert a Rotation to one with a smaller order! 
function Rotation{S1}(rot::Rotation{S2, N}) where {S1, S2, N}
    S1 > S2 && throw(DimensionMismatch(
        "Cannot convert a `Rotation` of order $S2 to order $S1"))
    Rotation(rot.m[1:S1])
end

function Rotation(m::DCM{N}, ω::AbstractVector) where N
    dm = ddcm(m, SVector(ω))
    return Rotation((m, dm))
end

# ---
# Type Conversions and Promotions 

# Static Arrays API 
Size(::Rotation{S, N}) where {S, N} = Size((3*S, 3*S))

similar_type(::Rotation{S, N}) where {S, N} = Rotation{S, N}
similar_type(::Rotation{S, <:Any}, ::Type{N}) where {S, N} = Rotation{S, N}

# Convert a Rotation to a tuple 
@generated function Base.Tuple(rot::Rotation{S, N}) where {S, N} 

    args = []
    for j = 1:3S 
        Oⱼ = (j-1) ÷ 3 + 1
        for i = 1:3S
            Oᵢ = (i-1) ÷ 3 + 1
            if Oⱼ > Oᵢ
                push!(args, N(0))
            else 
                row = i - 3*(Oᵢ - 1)
                col = j - 3*(Oⱼ - 1)

                push!(args, :(rot[$Oᵢ-$Oⱼ+1][$row, $col]))
            end
        end
    end

    return quote 
        Base.@_inline_meta
        @inbounds tuple($(args...))
    end
end

# Generic Rotation-to-StaticArrays conversions
@inline function (::Type{SA})(rot::Rotation) where {SA <: StaticArray}
    SA(Tuple(rot))
end

@inline MMatrix(rot::Rotation{S, N}) where {S, N} = MMatrix{3*S, 3*S}(rot)
@inline SMatrix(rot::Rotation{S, N}) where {S, N} = SMatrix{3*S, 3*S}(rot)

# Returns the inner datatype of a given DCM 
_dcm_type(::Union{DCM{T}, Type{DCM{T}}}) where T = T 

# Returns a promoted type for a given tuple of DCMs 
@generated function promote_dcm_eltype(::Union{T, Type{T}}) where T <: Tuple 
    t = Union{}
    for i = 1:length(T.parameters)
        tmp = _dcm_type(Base.unwrapva(T.parameters[i]))
        t = :(promote_type($t, $tmp))
    end 

    return quote 
        Base.@_inline_meta
        $t
    end
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
@inline Base.:*(A::Rotation{S, <:Any}, B::Rotation{S, <:Any}) where S = _compose_rot(A, B)

function Base.:*(::Rotation{S1, <:Any}, ::Rotation{S2, <:Any}) where {S1, S2}
    throw(DimensionMismatch("Cannot multiply two `Rotation` types of order $S1 and $S2"))
end

@generated function _compose_rot(A::Rotation{S, <:Any}, B::Rotation{S, <:Any}) where S
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
        DimensionMismatch("Cannot multiply `Rotation` of size ($sa, $sa) and a ($sb,) length vector"))

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
