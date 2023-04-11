module Utils

using ForwardDiff
using ForwardDiff: derivative
using LinearAlgebra
using ReferenceFrameRotations: DCM
using StaticArrays: SMatrix, SA, SVector


using InterfacesUtils
using InterfacesUtils: @filetype, AbstractFile

# IO 
include("IO/tpc.jl")

include("constants.jl")

include("geodesy.jl")

# Math
include("Math/derivatives.jl")
include("Math/vectors.jl")
include("Math/rotation.jl")

end
