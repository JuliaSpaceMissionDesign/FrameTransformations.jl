using ReferenceFrameRotations
using StaticArrays


struct Rotation{S, N}
    m::NTuple{S, DCM{N}}
end

# questa roba qua non piace a code_warntype e quindi fa l'allocco 
# function Rotation{N}(args...) where N
#     Rotation{length(args), N}(args)
# end

function Rotation(m::DCM{N}, ω::AbstractVector) where N
    dm = ddcm(m, SVector(ω))
    return Rotation{2, N}((m, dm))
end

function Rotation(m::DCM{N}) where N
    Rotation((m,))
end

function Rotation(m::DCM{N}, dm::DCM{N}) where N
    Rotation((m, dm))
end

