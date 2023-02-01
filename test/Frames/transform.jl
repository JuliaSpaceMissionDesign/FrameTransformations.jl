using SPICE: kclear, furnsh, sxform, bodc2n


@testset "Transformations" verbose=true begin 
    include("Definitions/planets.jl")
    include("Definitions/ecliptic.jl")
    include("Definitions/moon.jl")
    include("Definitions/earth.jl")
end;