using StaticArrays
using FrameTransformations

F = FrameSystem{2,Float64}()

add_axes!(F, :SatFrame, -1)

add_point!(F, :SC, -10000, :SatFrame)

sa_offset_left = [1.0, 0.0, 0.0]
sa_offset_right = [-1.0, 0.0, 0.0]
an_offset = [0.0, 0.0, -1.0]

add_point_fixedoffset!(F, :SolArrLeft, -10101, :SC, :SatFrame, sa_offset_left)
add_point_fixedoffset!(F, :SolArrRight, -10102, :SC, :SatFrame, sa_offset_right)
add_point_fixedoffset!(F, :Antenna, -10001, :SC, :SatFrame, an_offset)

vector3(F, :SolArrLeft, :SC, :SatFrame, 123.0)

vector6(F, :Antenna, :SolArrRight, :SatFrame, 456.0)

fun(t) = SA[cos(t), sin(t), 0.0]

add_point_dynamical!(F, :TimedAppendage, -10003, :SolArrLeft, :SatFrame, fun)

vector3(F, :TimedAppendage, :SC, :SatFrame, π / 3)

fun(t) = SA[cos(t), sin(t), 0]
dfun(t) = SA[cos(t), sin(t), 0, -sin(t), cos(t), 0]

add_point_dynamical!(F, :TimedAppendage2, -10004, :SolArrLeft, :SatFrame, fun, dfun)

vector6(F, :TimedAppendage2, :SC, :SatFrame, π / 3)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
