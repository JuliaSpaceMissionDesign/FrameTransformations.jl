"""
    AbstractFrameTransformation

An abstract type representing all frames transformation primitives.
"""
abstract type AbstractFrameTransformation end

"""
    AbstractFrameRotation{F1, F2}

An abstract type representing frames rotations
"""
abstract type AbstractFrameRotation{F1<:AbstractFrame, 
    F2<:AbstractFrame} <: AbstractFrameTransformation end

struct Rotation{F1,F2,T<:AbstractFloat} <: AbstractFrameRotation{F1, F2}
    origin::F1
    target::F2
    m::DCM{T}
    dm::DCM{T}
    function Rotation(
        origin::F1, 
        target::F2, 
        m::DCM{T},
        dm::DCM{T}=DCM(@SMatrix(zeros(T, 3, 3)))) where {F1, F2, T}
        new{F1, F2, T}(origin, target, m, dm)
    end
end

function Rotation(origin::F1, target::F2, m::DCM{T}, 
    ω::AbstractVector{T}) where {F1<:AbstractFrame, F2<:AbstractFrame, T<:AbstractFloat}
    dm = ddcm(m, SVector(ω))
    return Rotation(origin, target, m, dm)
end

function Rotation(origin::F, target::F, args...) where {F<:AbstractFloat}
    Rotation(origin, target, I3)
end

from(rot::Rotation) = rot.origin
to(rot::Rotation) = rot.target

function Base.inv(rot::Rotation)
    return Rotation(from(rot), to(rot), rot.m', rot.dm')
end

function compose(rot1::Rotation{F1, F}, rot2::Rotation{F, F2}) where {F<:AbstractFrame, 
    F1<:AbstractFrame, F2<:AbstractFrame}
    Rotation(from(rot1), to(rot2), rot2.m * rot1.m, rot2.dm * rot1.m + rot2.m * rot1.dm)
end

∘(rot1::Rotation, rot2::Rotation) = compose(rot1, rot2)

function apply(rot::Rotation, pos::VN, vel::VN) where {VN<:AbstractVector} 
    posn = rot.m * pos 
    veln = rot.dm * pos + rot.m * vel
    return posn, veln
end

function apply(rot::Rotation, vec::VN) where {VN<:AbstractVector} 
    n = length(vec)
    n < 3 && throw(ArgumentError("`vec` must have at least 3 elements."))
    pos = @view vec[1:3]
    posn = rot.m * pos
    n < 6 && return posn
    vel = @view vec[4:6]
    veln = rot.dm * pos + rot.m * vel
    return SA[posn[1], posn[2], posn[3], veln[1], veln[2], veln[3]]
end