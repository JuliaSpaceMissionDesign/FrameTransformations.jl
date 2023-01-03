module Orient

    using StaticArrays
    using LinearAlgebra
    using ReferenceFrameRotations
    using DelimitedFiles 
    using RemoteFiles 

    using Basic.Tempo 
    using Basic.Utils: arcsec2rad
    using Basic.Utils: InterpolationAkima, interpolate

    include("Earth/Earth.jl")
    include("Earth/eop.jl")
end