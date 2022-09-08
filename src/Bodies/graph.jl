export BodyGraph

const BODY_TYPES = (
    :CelestialBody,
    :Barycenter,
    :Planet,
    :NaturalSatellite,
    :MinorBody,
    :Asteroid,
    :Comet,
)

"""
    BodyGraph{T} where {T <: Integer}

Type alias of `NodeGraph{T, T} where {T <: Integer}`, to represent graphs of 
connected bodies.
"""
const BodyGraph{T} = NodeGraph{T, T} where {T <: Integer}

"""
    BodyGraph()

Convenience constructor for graphs of bodies. 
Return a BodyGraph of `NAIFId` body types.
"""
function BodyGraph(::Type{T}) where {T <: Integer}
    BodyGraph{T}(SimpleGraph())
end
function BodyGraph()
    BodyGraph(Int64)
end

"""
    register!(g::BodyGraph{N}, id::N) where {N <: Integer}

Register a new object with identifier `id` in the body graph `g`.
"""
register!(g::BodyGraph{N}, id::N) where {N <: Integer} = add_vertex!(g, id)

"""
    connect!(g::BodyGraph{N}, id1::N, id2::N) where {N <: Integer}

Connect two object within a body graph `g`.
"""
connect!(g::BodyGraph{N}, id1::N, id2::N) where {N <: Integer} = add_edge!(g, id1, id2)

"""
    find_path(g::BodyGraph, from::CelestialBody, to::CelestialBody)
    find_path(g::BodyGraph, from::N, to::N) where {N <: Integer}

Find the shortest path linking two objects in a BodyGraph-type graph.
"""
function find_path(g::BodyGraph{N}, from::CelestialBody, 
    to::CelestialBody) where {N <: Integer}
    find_path(g, body_naifid(from), body_naifid(to))
end

function find_path(g::BodyGraph, from::N, to:: N) where {N <: Integer}
    get_nodes(g, from, to)
end