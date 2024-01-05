using FrameTransformations

G = FrameSystem{4, Float64}()

add_axes_inertial!(G, :ICRF, 1)
add_point_root!(G, :Earth, 399, 1)

using JSMDInterfaces.Math: interpolate
using JSMDUtils.Math: InterpCubicSplines

c = 2π/86400
de = 0:3600:2*86400;
fv(t) = [cos(c*t), sin(c*t), 0, -c*sin(c*t), c*cos(c*t), 0];
fv2(t) = [sin(c*t), cos(c*t), 0];
const ORB1 = InterpCubicSplines(de, hcat([fv(e) for e in de]...));
const ORB2 = InterpCubicSplines(de, hcat([fv2(e) for e in de]...));

add_point_dynamical!(G, :SC1, -1, 399, 1, t->interpolate(ORB1, t)[1:3], t->interpolate(ORB1, t))
add_point_dynamical!(G, :SC2, -2, 399, 1, t->interpolate(ORB2, t))

vector12(G, 399, -1, 1, 12345.0)

vector12(G, -1, -2, 1, 12345.0)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
