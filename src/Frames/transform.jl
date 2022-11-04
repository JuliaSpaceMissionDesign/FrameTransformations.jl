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
    
"""
    Rotation{F1, F2, T<:AbstractFloat}

A type representing Rotations between reference frames.
"""
struct Rotation{F1, F2, T<:AbstractFloat} <: AbstractFrameRotation{F1, F2}
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

frame_origin(rot::Rotation) = rot.origin
frame_target(rot::Rotation) = rot.target

function Base.inv(rot::Rotation)
    return Rotation(frame_origin(rot), frame_target(rot), rot.m', rot.dm')
end

function compose(rot1::Rotation{F1, F}, rot2::Rotation{F, F2}) where {F<:AbstractFrame, 
    F1<:AbstractFrame, F2<:AbstractFrame}
    Rotation(
        frame_origin(rot1), 
        frame_target(rot2), 
        rot2.m * rot1.m, rot2.dm * rot1.m + rot2.m * rot1.dm
    )
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

function Rotation(from::F1, to::F2, ep::Epoch) where {F1<:AbstractFrame,F2<:AbstractFrame}
    frames = find_path(from, to)
    return _Rotation(ep, frames...)
end

function _Rotation(ep::Epoch, path::AbstractFrame...)::Rotation
    f1 = path[1]
    f2 = path[2]
    rot = Rotation(f1, f2, ep)
    for i in 2:length(path)-1
        f1 = path[i]
        f2 = path[i+1]
        rot = compose(rot, Rotation(f1, f2, ep))
    end
    return rot
end