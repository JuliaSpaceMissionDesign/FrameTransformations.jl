using ReferenceFrameRotations
using StaticArrays


struct Rotation{S, N}
    m::NTuple{S, DCM{N}}
end

Rotation(m::DCM) = Rotation((m,))

function Rotation(m::DCM{N}, dm::DCM{N}) where N
    Rotation((m, dm))
end

function Rotation(m::DCM{N}, dm::DCM{N}, ddm::DCM{N}) where N 
    Rotation((m, dm, ddm))
end

function Rotation(m::DCM{N}, ω::AbstractVector) where N
    dm = ddcm(m, SVector(ω))
    return Rotation{2, N}((m, dm))
end

