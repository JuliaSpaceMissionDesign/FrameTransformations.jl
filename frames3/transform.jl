using ReferenceFrameRotations
using StaticArrays
using ForwardDiff
using BenchmarkTools

import LinearAlgebra.matprod

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

# Constructor with angular velocity 
function Rotation(m::DCM{N}, ω::AbstractVector) where N
    dm = ddcm(m, SVector(ω))
    return Rotation{2, N}((m, dm))
end

Base.getindex(R::Rotation, i) = R.m[i]

@inline Base.inv(rot::Rotation) = _inverse_rot(rot)
@generated function _inverse_rot(rot::Rotation{S, N}) where {S, N}
    return quote 
        @inbounds Rotation($((:(rot.m[$i]') for i in 1:S)...))
    end
end


@inline Base.:*(A::Rotation, b::AbstractVector) = _mul_rot(A, b)
@inline Base.:*(A::Rotation{S, <:Any}, B::Rotation{S, <:Any}) where S = _multiply_rot(A, B)

function Base.:*(::Rotation{S1, <:Any}, ::Rotation{S2, <:Any}) where {S1, S2}
    throw(ArgumentError("Cannot multiply two `Rotation` of different order!"))
end


@generated function _multiply_rot(A::Rotation{S, N}, B::Rotation{S, N}) where {S, N}

    expr = Expr(:call, :Rotation)
    
    for i in 1:S
        sum_expr = Expr(:call, :+)
        for j in 1:i 
            c = binomial(i-1, j-1)
            ai = Expr(:ref, :A, i-j+1)
            bi = Expr(:ref, :B, j)

            push!(sum_expr.args, Expr(:call, :*, c, ai, bi))
        end
        push!(expr.args, sum_expr)    
    end

    return quote 
        @inbounds $(expr)
    end
end

# Function to compute product between Rotation and a generic vector
@generated function _mul_rot(A::Rotation{S, Ta}, b::AbstractVector{Tb}) where {S, Ta, Tb}

    exprs = [[Expr(:call, :+) for _ = 1:3] for _ = 1:S]

    for i = 1:S 
        for j = 1:i
            for k = 1:3 
                push!(exprs[i][k].args, 
                    StaticArrays.combine_products([:(A[$i-$j+1][$k, $w]*b[3*($j-1)+$w]) for w = 1:3]))
            end 
        end
    end

    sa = 3*S 
    retexpr = :(@inbounds return similar_type(b, T, Size($sa))(tuple($((exprs...)...))))
    return quote 
        Base.@_inline_meta
        length(b) != $sa && throw(DimensionMismatch(
            "Tried to multiply arrays of size ($($sa), $($sa)) and $(size(b))"))

        T = Base.promote_op(matprod, Ta, Tb)
        $retexpr

    end
end

# Function to compute product between Rotation and a static vector
@generated function _mul_rot(A::Rotation{S, Ta}, b::StaticVector{sb, Tb}) where {S, sb, Ta, Tb}
    sa = 3*S 
    sa != sb && throw(DimensionMismatch("Tried to multiply arrays of size ($sa, $sa) and ($sb,)"))

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
        T = Base.promote_op(matprod, Ta, Tb)
        $retexpr
    end
end

