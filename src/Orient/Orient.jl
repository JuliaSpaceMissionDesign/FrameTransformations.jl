module Orient

    using StaticArrays
    using LinearAlgebra
    using ReferenceFrameRotations

    using Basic.Tempo 
    using Basic.Utils: arcsec2rad

    include("Earth/Earth.jl")
    include("Earth/eopdata.jl")
end