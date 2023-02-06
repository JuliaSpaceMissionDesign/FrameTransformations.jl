using SPICE: kclear, furnsh, sxform, bodc2n


@testset "Transformations" verbose=true begin 
    
    include("axes.jl")
    include("points.jl")
    
    # FIXME: test planets da ripristinare
    # include("Definitions/planets.jl")
    include("Definitions/ecliptic.jl")
    include("Definitions/moon.jl")
    include("Definitions/earth.jl")

end;