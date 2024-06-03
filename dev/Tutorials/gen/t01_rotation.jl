using FrameTransformations
using ReferenceFrameRotations

dcm  = angle_to_dcm(π/3, :Z)
δdcm = DCM(0I)

R = Rotation(dcm, δdcm)

R[1]

R[2]

v = [1., -6., 3., 0., 5., 0]
R*v

inv(R)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
