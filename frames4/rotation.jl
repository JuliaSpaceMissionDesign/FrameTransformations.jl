using ReferenceFrameRotations
using StaticArrays
using LinearAlgebra: matprod, UniformScaling

import StaticArrays: similar_type, Size

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

Base.size(::Rotation{S, <:Any}) where S = (3S, 3S)
Base.Tuple(rot::Rotation) = rot.m
Base.getindex(R::Rotation, i) = R.m[i]

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

function Rotation(m::DCM{N}, ω::AbstractVector) where N
    dm = ddcm(m, SVector(ω))
    return Rotation((m, dm))
end

# ---
# Type Operations and Promotions 

# Static Arrays API 
Size(::Rotation{S, N}) where {S, N} = Size((3*S, 3*S))

similar_type(::Rotation{S, N}) where {S, N} = Rotation{S, N}
similar_type(::Rotation{S, <:Any}, ::Type{N}) where {S, N} = Rotation{S, N}

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
        Base.Base.@_inline_meta
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
    throw(
        ArgumentError("Cannot multiply two `Rotation` of different order!")
    )
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
