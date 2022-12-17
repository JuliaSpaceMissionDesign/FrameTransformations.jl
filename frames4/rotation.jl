using ReferenceFrameRotations
using StaticArrays
using ForwardDiff
using BenchmarkTools

struct Rotation{S, N}
    m::NTuple{S, DCM{N}}
end

Rotation(m::DCM{N}) where N = Rotation((m,))
Rotation(m::DCM{N}, dm::DCM{N}) where N = Rotation((m, dm))
Rotation(m::DCM{N}, dm::DCM{N}, ddm::DCM{N}) where N = Rotation((m, dm, ddm))

function Rotation(m::DCM{N}, ω::AbstractVector) where N
    dm = ddcm(m, SVector(ω))
    return Rotation{2, N}((m, dm))
end

Base.inv(rot::Rotation) = _inverse_rot(rot)
@generated function _inverse_rot(rot::Rotation{S, N}) where {S, N}
    quote 
        @inbounds Rotation($((:(rot.m[$i]') for i in 1:S)...))
    end
end

Base.:*(A::Rotation{S, N}, B::Rotation{S, N}) where {S, N} = _multiply_rot(A, B)
@generated function _multiply_rot(A::Rotation{S, N}, B::Rotation{S, N}) where {S, N}
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