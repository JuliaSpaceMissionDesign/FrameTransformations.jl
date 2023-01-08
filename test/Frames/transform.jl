using SPICE: kclear, furnsh, sxform, bodc2n


@testset "Transformations" verbose=true begin 
    include("Definitions/planets.jl")
    include("Definitions/ecliptic.jl")
end;