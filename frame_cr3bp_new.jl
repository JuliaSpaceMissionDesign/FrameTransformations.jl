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

function dtest(node32, t, x, y)
    node32.dfun(t, x, y)
end

function ddtest(node, t, x, y)
    node.ddfun(t, x, y)
end

ta = 2.
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

FRAMES = FrameSystem{Float64}();
FRAMES2 = FrameSystem{Float64}();

add_inertial_axes!(FRAMES, :IRF, 1)
add_inertial_axes!(FRAMES2, :IRF, 1)
add_rotating_axes!(FRAMES, :Synodic, 5, IRF, in2syn)
add_rotating_axes!(FRAMES2, :Synodic, 5, IRF, in2syn, din2syn)

in2syn(t::T) where T = angle_to_dcm(-t, :Z)

node = FRAMES.axes_graph.nodes[end]
node2 = FRAMES2.axes_graph.nodes[end]

@benchmark test($node, $ta, $xa, $ya)
@benchmark test($node2, $ta, $xa, $ya)

@benchmark dtest($node, $ta, $xa, $ya)
@benchmark dtest($node2, $ta, $xa, $ya)

# @benchmark ddtest($node, $t, $x, $y)

# @code_warntype dtest(node, t, x, y)