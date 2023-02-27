using SPICE: kclear, furnsh, sxform, bodc2n
using Basic.Utils: D¹, D², D³, angle_to_δdcm, angle_to_δ²dcm, angle_to_δ³dcm


@testset "Transformations" verbose=true begin 
    
    include("axes.jl")
    include("points.jl")
    
    include("Definitions/planets.jl")
    include("Definitions/ecliptic.jl")
    include("Definitions/moon.jl")

end;