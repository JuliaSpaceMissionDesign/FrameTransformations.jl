export BodyGraph

import Basic: register!, connect!, find_path

function BodyGraph(::Type{N}) where {N<:Integer}
    NodeGraph{NAIFId, N, N}(SimpleGraph{N}())
end

BodyGraph() = BodyGraph(Int64)

"""
    register!(g::NodeGraph, id::Integer) 
    register!(g::NodeGraph, id::NAIFId)

Register a new object with identifier `id` in the bodies graph.
"""
register!(g::NodeGraph, id::NAIFId) = add_vertex!(g, id)
register!(g::NodeGraph, id::Integer) = register!(g, NAIFId(id))

"""
    connect!(g::NodeGraph, id1::Integer, id2::Integer)
    connect!(g::NodeGraph, id1::NAIFId, id2::NAIFId)

Connect two object within the bodies graph.
"""
connect!(g::NodeGraph, id1::NAIFId, id2::NAIFId) = add_edge!(g, id1, id2)
connect!(g::NodeGraph, id1::Integer, id2::Integer) = connect!(g, NAIFId(id1), NAIFId(id2))

"""
    find_path(g::NodeGraph, from::CelestialBody, to::CelestialBody)
    find_path(g::NodeGraph, from::N, to::N) where {N <: Integer}

Find the shortest path linking two objects in a `NodeGraph`-type graph.
"""
function find_path(g::NodeGraph, from::CelestialBody, to::CelestialBody) 
    find_path(g, body_naifid(from), body_naifid(to))
end

function find_path(g::NodeGraph, from::NAIFId, to::NAIFId)
    get_nodes(g, from, to)
end

