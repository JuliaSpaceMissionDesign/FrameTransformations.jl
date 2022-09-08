export BODIES 

"""
    register!(id::Integer) 
    register!(id::NAIFId)

Register a new object with identifier `id` in the bodies graph.
"""
register!(id::NAIFId) = add_vertex!(BODIES, id)
register!(id::Integer) = register!(NAIFId(id))

"""
    connect!(id1::Integer, id2::Integer)
    connect!(id1::NAIFId, id2::NAIFId)

Connect two object within the bodies graph.
"""
connect!(id1::NAIFId, id2::NAIFId) = add_edge!(BODIES, id1, id2)
connect!(id1::Integer, id2::Integer) = connect!(id1, id2)

"""
    find_path(g::BodyGraph, from::CelestialBody, to::CelestialBody)
    find_path(g::BodyGraph, from::N, to::N) where {N <: Integer}

Find the shortest path linking two objects in a BodyGraph-type graph.
"""
function find_path(from::CelestialBody, to::CelestialBody) 
    find_path(BODIES, body_naifid(from), body_naifid(to))
end

function find_path(from::NAIFId, to::NAIFId)
    get_nodes(BODIES, from, to)
end

const BODIES = NodeGraph{CelestialBody, NAIFId, Int64}(SimpleGraph{Int64}())