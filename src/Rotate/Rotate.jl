module Rotate

import Base: +, -, *, /, \, ≈, ==, ∘
import Base: conj, convert, copy, display, eltype, firstindex, getindex, imag
import Base: inv, iterate, lastindex, ndims, one, rand, real, setindex!, show
import Base: summary, zero, zeros
import Base: Broadcast.broadcastable
import LinearAlgebra: norm
import StaticArrays: similar_type

using LinearAlgebra
using StaticArrays

# Re-export `I` from LinearAlgebra.
export I

# TODO: review and test API

include("types.jl")

include("dcm.jl")
include("angle.jl")
include("angleaxis.jl")
include("compose_rotations.jl")
include("inv_rotations.jl")
include("quaternion.jl")

include("conversions/angle_to_angle.jl")
include("conversions/angle_to_angleaxis.jl")
include("conversions/angle_to_dcm.jl")
include("conversions/angle_to_quat.jl")
include("conversions/angle_to_rot.jl")
include("conversions/angleaxis_to_angle.jl")
include("conversions/angleaxis_to_dcm.jl")
include("conversions/angleaxis_to_quat.jl")
include("conversions/api.jl")
include("conversions/dcm_to_angle.jl")
include("conversions/dcm_to_angleaxis.jl")
include("conversions/dcm_to_quat.jl")
include("conversions/quat_to_angle.jl")
include("conversions/quat_to_angleaxis.jl")
include("conversions/quat_to_dcm.jl")
include("conversions/smallangle_to_dcm.jl")
include("conversions/smallangle_to_quat.jl")
include("conversions/smallangle_to_rot.jl")

end