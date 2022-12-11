include("../src/errors.jl")

import Graphs: 
    AbstractGraph, 
    SimpleGraph, 
    add_edge!,
    add_vertex!, 
    nv, 
    dijkstra_shortest_paths, 
    enumerate_paths, 
    has_vertex, 
    has_path    

abstract type AbstractGraphNode end 

# Function to return id associated to given node 
function get_node_id(b::AbstractGraphNode)
    throw(NotImplementedError(
        "`get_id` shall be implemented for $b")
        )
end

struct GraphSystem{N <: AbstractGraphNode, G}
    graph::G             
    ids::Dict{Int, Int}                      # mapping between IDs (NAIF) and nodes idx 
    nodes::Vector{N}                         # nodes 
    paths::Dict{Int, Dict{Int, Vector{Int}}} # mapping between NAIF (NAIF and nodes idx)
    edges::Dict{Int, Dict{Int, Int}}         # edges  

    function GraphSystem{N}(g::G) where {G <: AbstractGraph, N}
        new{N, G}(g, 
            Dict{Int, Int}(), 
            Vector{N}(undef, 0),
            Dict{Int, Dict{Int, Int}}(), 
            Dict{Int, Dict{Int, Int}}(),
        )
    end
end
GraphSystem{N}(g=SimpleGraph()) where {N} = GraphSystem{N}(g)

function GraphSystem(::Type{N}) where {N}
    GraphSystem{N}(SimpleGraph{Int}())
end

# get_iid --> ritorna l'internal ID del nodo avendo come input NAIFID 
# get_node --> restituisce axis con input internal ID 

get_iid(g::GraphSystem, node::Int) = g.ids[node] 
@inbounds get_node(g::GraphSystem, iid::Int) = g.nodes[iid] 
get_node_from_id(g::GraphSystem, node::Int) = get_node(g, get_iid(g, node))

# true if node is within graph g 
has_vertex(g::GraphSystem, node::Int) = haskey(g.ids, node)

# true if there is a path between `from` and `to`
has_path(g::GraphSystem, from::Int, to::Int) = has_path(g.graph, get_iid(g, from), get_iid(g, to))

Base.isempty(g::GraphSystem) = Base.isempty(g.nodes)

function add_vertex!(g::GraphSystem{N}, node::N) where {N <: AbstractGraphNode}
    id = get_node_id(node)
    has_vertex(g, id) && return 
    
    add_vertex!(g.graph)

    # Internal ID value 
    iid = nv(g.graph) 

    # Updates internal vectors\dictionaries
    push!(g.ids, id => iid)
    push!(g.nodes, node)
    nothing
end


function add_edge!(g::GraphSystem{N}, from::Int, to::Int, 
                   cost::Int=0) where {N}

    # Ensure the two vertexes already exist in the graph 
    if !(has_vertex(g, from) && has_vertex(g, to))
        throw(ErrorException("The vertex $from or $to is not contained in the graph."))
    end

    add_edge!(g.graph, get_iid(g, from), get_iid(g, to))
    compute_paths!(g)

    edges = get!(g.edges, from, Dict{Int, Int}())
    push!(edges, to => cost)
    nothing 
end

function get_edgecost(g::GraphSystem{N}, from::Int, to::Int) where {N}
    path = get_path(g, from, to)
    isempty(path) && return Int[]
    edges = Vector{Int}(undef, length(path)-1)
    for i in eachindex(edges)
        edges[i] = g.edges[path[i]][path[i+1]]
    end
    edges 
end

function get_path(g::GraphSystem{N}, from::Int, to::Int) where N
    (has_vertex(g, from) && has_vertex(g, to)) || return Int[]
    g.paths[from][to]
end

function compute_paths!(g::GraphSystem{N}) where {N}
    for (oiid, origin) in enumerate(g.nodes) 
        oid = get_node_id(origin)

        ds = dijkstra_shortest_paths(g.graph, oiid)
        for (tiid, target) in enumerate(g.nodes)
            oiid == tiid && continue 
            tid = get_node_id(target)

            path = enumerate_paths(ds, tiid)
            paths = get!(g.paths, oid, Dict{Int, Vector{Int}}())
            push!(paths, tid => path)
        end
    end

    nothing 
end

function connect!(g::GraphSystem, from::Int, to::Int)
    add_edge!(g, from, to)
end


