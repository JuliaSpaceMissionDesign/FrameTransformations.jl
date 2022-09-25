export BodyGraph

import Basic: register!, connect!, find_path

"""
    BodyGraph 

Convenience constructor for `NodeGraph`
"""
function BodyGraph(::Type{N}, ::Type{B}) where {N<:Integer, B}
    NodeGraph{B, N, N}(SimpleGraph{N}())
end
BodyGraph(::Type{B}) where {B} = BodyGraph(Int64, B)
BodyGraph() = BodyGraph(NAIFId)

"""
    register!(g::NodeGraph, id::Integer) 
    register!(g::NodeGraph, id::NAIFId)

Register a new object with identifier `id` in the bodies graph.
"""
function register!(g::NodeGraph{NAIFId, N, G, N}, id::NAIFId) where {N<:Integer, G} 
    add_vertex!(g, id)
end
function register!(g::NodeGraph{NAIFId, N, G, N}, id::Integer) where {N<:Integer, G}
    register!(g, NAIFId(id))
end
function register!(g::NodeGraph{B, N, G, N}, b) where {B<:CelestialBody, N<:Integer, G} 
    add_vertex!(g, b)
end

"""
    connect!(g::NodeGraph, id1::Integer, id2::Integer)
    connect!(g::NodeGraph, id1::NAIFId, id2::NAIFId)

Connect two object within the bodies graph.
"""
function connect!(g::NodeGraph{NAIFId, N, G, N}, 
    id1::NAIFId, id2::NAIFId) where {N<:Integer, G}
    add_edge!(g, id1, id2)
end
function connect!(g::NodeGraph{NAIFId, N, G, N}, 
    id1::Integer, id2::Integer) where {N<:Integer, G}
    connect!(g, NAIFId(id1), NAIFId(id2))
end
function connect!(g::NodeGraph{B, N, G, N}, 
    id1, id2) where {B<:CelestialBody, N<:Integer, G}
    connect!(g, id1, id2)
end

"""
    find_path(g::NodeGraph, from::CelestialBody, to::CelestialBody)
    find_path(g::NodeGraph, from::NAIFId, to::NAIFId)
    find_path(g::NodeGraph, from::N, to::N) where {N <: Integer}

Find the shortest path linking two objects in a `NodeGraph`-type graph.
"""
function find_path(g::NodeGraph{NAIFId, N, G, N}, 
    from::CelestialBody, to::CelestialBody) where {N<:Integer, G}
    find_path(g, body_naifid(from), body_naifid(to))
end

function find_path(g::NodeGraph{NAIFId, N, G, N}, 
    from::NAIFId, to::NAIFId) where {N<:Integer, G}
    get_nodes(g, from, to)
end

function find_path(g::NodeGraph{NAIFId, N, G, N}, 
    from::N, to::N) where {N<:Integer, G}
    find_path(g, NAIFId(from), NAIFId(to))
end

function find_path(g::NodeGraph{B, N, G, N}, 
    from, to) where {B<:CelestialBody, N<:Integer, G}
    find_path(g, from, to)
end
