module Bodies 

    import NodeGraphs: NodeGraph, SimpleDiGraph, SimpleGraph, add_edge!, add_vertex!, get_edgenodes
    
    export NAIFId, from_naifid, body_gm, body_mean_radius, body_equatorial_radius, body_polar_radius,
        body_parent
        
    include("types.jl")
    include("graph.jl")

end