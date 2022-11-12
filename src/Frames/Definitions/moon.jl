using CALCEPH: Ephem as CalcephEphemHandler, 
               unsafe_orient!, 
               unitRad, unitSec, useNaifId

using Basic.Utils: arcsec2rad

struct MoonPrincipalAxisFrame <: AbstractBodyCentricRotatingFrame 
    ptr::CalcephEphemHandler
    angles::Vector
    function MoonPrincipalAxisFrame(ptr::CalcephEphemHandler) 
        new(ptr, zeros(6))
    end
end

function build(::Type{MoonPrincipalAxisFrame}, ptr::CalcephEphemHandler, args...)
    return @eval begin
        const MOON_PA = MoonPrincipalAxisFrame($ptr)
        export MOON_PA
    end
end

function Rotation(
    origin::InternationalCelestialReferenceFrame, target::MoonPrincipalAxisFrame,
    e::Epoch 
)
    unsafe_orient!(
        target.angles, target.ptr, 
        DJ2000, j2000(e), 
        31008, # Moon Principal Axis NAIF ID
        unitRad+unitSec+useNaifId, 1
    )
    ϕ, θ, ψ = @view target.angles[1:3]
    R = angle_to_dcm(ϕ, θ, ψ, :ZXZ)
    ω = SA[target.angles[4], target.angles[5], target.angles[6]]
    Rotation(origin, target, R, ω)
end

function Rotation(
    origin::MoonPrincipalAxisFrame, target::InternationalCelestialReferenceFrame,
    e::Epoch
)
    inv(Rotation(target, origin, e))
end

struct MoonMeanEquatorFrame <: AbstractBodyCentricRotatingFrame end 

const MOON_ME = MoonMeanEquatorFrame()

# Rotation martix defined as for DE440
# https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/Tutorials/pdf/individual_docs/23_lunar-earth_pck-fk.pdf
# https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_de440_220930.tf
const _ROTMAT_PA2ME = angle_to_dcm(
    arcsec2rad.((-67.8526, -78.6944, -0.2785))..., 
    :ZYX
)

function Rotation(origin::MoonPrincipalAxisFrame, target::MoonMeanEquatorFrame, e::Epoch)
    Rotation(origin, target, _ROTMAT_PA2ME)
end

function Rotation(origin::MoonMeanEquatorFrame, target::MoonPrincipalAxisFrame, e::Epoch)
    Rotation(origin, target, _ROTMAT_PA2ME')
end