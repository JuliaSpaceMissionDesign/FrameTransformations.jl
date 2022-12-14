using BenchmarkTools 
using ForwardDiff

# CR3BP scheme 
include("frames/cph.jl")
include("frames/graph.jl")
include("frames3/transform.jl")
include("frames3/types.jl")
include("frames3/axes.jl")

function test(node, t, x, y)
    node.fun(t, x, y)
end

function dtest(node, t, x, y)
    node.dfun(t, x, y)
end

function ddtest(node, t, x, y)
    node.ddfun(t, x, y)
end

ta = 2.0
xa = SA[1., 2., 3.]
ya = SA[2., 3., 4.]

@axes IRF 1 InertialReferenceFrame
@axes MEME2000 2 
@axes ECLIPJ2000 4 EclipticJ2000

in2syn(t::T) where T = angle_to_dcm(-t, :Z)

function din2syn(t::T) where T 
    s, c = sincos(t)
    angle_to_dcm(-t, :Z), DCM{T}(-s, -c, 0, c, -s, 0, 0, 0, 0)
end

in2syn2(t::T, x, y) where T = Rotation(angle_to_dcm(-t, :Z))

FRAMES = FrameSystem{Float64}();

add_inertial_axes!(FRAMES, :IRF, 1)
add_rotating_axes!(FRAMES, :Synodic, 5, IRF, in2syn, din2syn)

in2syn(t::T) where T = angle_to_dcm(-t, :Z)

node1 = FRAMES.axes_graph.nodes[1]
node2 = FRAMES.axes_graph.nodes[end]

@benchmark test($node1, $ta, $xa, $ya)
@benchmark test($node2, $ta, $xa, $ya)

# @benchmark dtest($node, $ta, $xa, $ya)
# @benchmark ddtest($node, $ta, $xa, $ya)

# @benchmark ddtest($node, $t, $x, $y)

# @code_warntype dtest(node, t, x, y)