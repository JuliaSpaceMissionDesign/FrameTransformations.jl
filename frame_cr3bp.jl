using BenchmarkTools 

# CR3BP scheme 
include("frames/Frames.jl")

FRAMES = FrameSystem{Float64}();

in2syn(t::T, y) where T = angle_to_dcm(-t, :Z)

function din2syn(t::T, y) where T 
    s, c = sincos(t)
    DCM{T}(-s, -c, 0, c, -s, 0, 0, 0, 0)
end

@inertial_axes FRAMES IRF 1 type=InertialReferenceFrame
@rotating_axes FRAMES SRF 2 IRF in2syn din2syn type=SynodicReferenceFrame

μ = 0.0153 

@root_point FRAMES emb 3 SRF type=EarthMoonBarycenterPoint
@fixed_point FRAMES moon 301 emb SRF SA[1-μ, 0., 0.] type=MoonPoint
@fixed_point FRAMES earth 399 emb SRF SA[-μ, 0., 0.] type=EarthPoint

@benchmark Rotation($FRAMES, $IRF, $SRF, $(0.1))


get_vector3(FRAMES, moon, earth, IRF, 0.1)

@benchmark get_vector3($FRAMES, $earth, $moon, $IRF, $(0.1))