using ReferenceFrameRotations
import FunctionWrappers: FunctionWrapper
using StaticArrays
using ForwardDiff

## Types

abstract type AbstractGeometricalTransformation{T} end

struct Rotation{T<:AbstractFloat, S} <: AbstractGeometricalTransformation{T}
    m::NTuple{S, DCM{T}}
end

struct Translation{T<:AbstractFloat, S} <: AbstractGeometricalTransformation{T}
    v::MVector{S, T}
end


# Points 

struct FramePoint{T, N}
    V::Vector{Translation{T, N}}
    fin::FunctionWrapper{Nothing, Tuple{Translation{T, N}, T}}
end


# Axes

struct FrameAxes{T, N, M}
    R::Vector{Rotation{T, N}}
    fop::FunctionWrapper{Rotation{T, N}, Tuple{T, SVector{M, T}}}
end

# 0 pos, 1 pos vel, 2 pos vel acc
function FrameAxes(::Type{T}, order::N, fun) where {T<:AbstractFloat, N}
    return FrameAxes{T, order+1, 3*(order+1)}(
        Vector{Rotation{T, order+1}}(), fun
    ) 
end

function FrameAxes(order::N, fun) where {N}
    return FrameAxes(Float64, order, fun)
end
























function rot(t::T, x::SVector{3, T})::Rotation{T, 1} where {T<:AbstractFloat}
    return Rotation{T, 1}((DCM{T}(1.0I),))
end 

function rot(t::T, x::SVector{6, T})::Rotation{T, 2} where {T<:AbstractFloat}
    return Rotation{T, 1}((DCM{T}(1.0I), DCM{T}(1.0I)))
end 

a = FrameAxes(0, rot)

function compute(a::FrameAxes)
    return a.fun(0.0, SVector{3, Float64}(1.0, 2.0,3))
end