using BenchmarkTools
using ForwardDiff
using Basic 

# CR3BP scheme 
include("frames/Frames.jl")

eph = CalcephProvider(["/home/michele/spice/kernels/spk/de440.bsp"])

F = FrameSystem{Float64}(eph);


icrf2meme = angle_to_dcm(pi*2, :X)
meme2eclip = angle_to_dcm(pi*2, pi/3, :XY)

# NEW frames set-up
@axes ICRF 1 InternationalCelestialReferenceFrame
@axes MEME2000 2 MeanEarthMeanEquinoxJ2000
@axes ECLIPJ2000 3 EclipticEquinoxJ2000
@axes IAU_EARTH 4 IauEarth

add_axes_inertial!(F1, ICRF)
add_axes_inertial!(F1, MEME2000; parent=1, dcm=icrf2meme)
add_axes_inertial!(F1, ECLIPJ2000; parent=MEME2000, dcm=meme2eclip)
add_axes_inertial!(F1, IAU_EARTH; parent=MEME2000, dcm=meme2eclip)

using Logging
function test_log(x)
    @debug "debugger"
    @info "diologger"
    @warn "mado"
    @error "stoppa"
end

# OLD frames set-up 
# @inertial_axes FA  ICRFA       1 type=InternationalCelestialReferenceFrameA
# @inertial_axes FA  MEME2000A   2 type=MeanEquatorMeanEquinoxJ2000A parent=ICRFA dcm=icrf2meme
# @inertial_axes FA  ECLIPJ2000A 3 type=EclipticEquinoxJ2000A parent=MEME2000A dcm=meme2eclip 
# @inertial_axes FA  IAU_EARTHA  4 type=IauEarthA parent=ECLIPJ2000A dcm=meme2eclip 

# # new test with integer classes 
# add_axes_inertial!(FC, ICRF)
# add_axes_inertial!(FC, MEME2000; parent=1, dcm=icrf2meme)
# add_axes_inertial!(FC, ECLIPJ2000; parent=MEME2000, dcm=meme2eclip)
# add_axes_inertial!(FC, IAU_EARTH; parent=ECLIPJ2000, dcm=meme2eclip)

# # NEW vs OLD performance comparisons 
# @benchmark get_rotation3($F1, $ICRF, $MEME2000, 0.)
# @benchmark RotationA($FA, $ICRFA, $MEME2000A, 0.)

# @benchmark get_rotation3($F1, $MEME2000, $ECLIPJ2000, 0.)
# @benchmark RotationA($FA, $MEME2000A, $ECLIPJ2000A, 0.)

# # why is this slow?
# @benchmark get_rotation3($F1, $ICRF, $ECLIPJ2000, 0.)
# @benchmark RotationA($FA, $ICRFA, $ECLIPJ2000A, 0.)

n1 = F1.axes.nodes[1] 
n2 = F1.axes.nodes[2]

nA1 = FA.axes_graph.nodes[1]
nA2 = FA.axes_graph.nodes[2]

nC1 = FC.axes.nodes[1] 
nC2 = FC.axes.nodes[2]

path2 = [1, 2, 3]

@benchmark _RotationA($FA, 0., $path2)
@benchmark _compute_rot3($F1, 0., $path2)
@benchmark _compute_rot3($FC, 0., $path2)

@point Earth 399 EarthPoint 

@orient_iau_angles Earth 


struct EarthPoint 
    NAIFId::Val{399}
end


struct MarsPoint 
    NAIFId::Val{299}
end
