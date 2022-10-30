module Orient

    using StaticArrays
    using LinearAlgebra
    using ReferenceFrameRotations

    using Basic.Tempo 
    using Basic.Utils: arcsec2rad

    include("Planets/abstract.jl")
    include("Planets/angles.jl")

    include("Earth/types.jl")
    include("Earth/obliquity.jl")
    include("Earth/precession.jl")
    include("Earth/nutation.jl")
end