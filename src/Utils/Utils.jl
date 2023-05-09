module Utils

using LinearAlgebra
using ReferenceFrameRotations: DCM
using StaticArrays: SMatrix, SA, SVector

using SMDInterfacesUtils
import SMDInterfacesUtils.IO: load
using SMDInterfacesUtils.IO: @filetype, AbstractFile

# IO 
include("IO/tpc.jl")

include("constants.jl")
include("geodesy.jl")

# Math
include("Math/vectors.jl")
include("Math/rotation.jl")

end
