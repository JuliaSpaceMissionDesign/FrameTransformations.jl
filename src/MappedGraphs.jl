module MappedGraphs

    export AbstractGraphNode, 
           MappedNodeGraph, MappedGraph, MappedDiGraph,
           SimpleGraph, SimpleDiGraph,
           get_mappedid, get_mappednode, get_node_id, get_path,
           has_vertex, has_path, add_edge!, add_vertex!, get_edgecosts

    import Graphs: 
        AbstractGraph, 
        SimpleGraph, 
        SimpleDiGraph,

        add_edge!,
        add_vertex!, 
        nv, 
        dijkstra_shortest_paths, 
        enumerate_paths, 
        has_vertex, 
        has_path   


    abstract type AbstractGraphNode end 

    function get_node_id(b::AbstractGraphNode)
        throw(ErrorException("`get_node_id` shall be implemented for $b"))
    end

    struct MappedNodeGraph{N, G}
        graph::G             
        mid::Dict{Int, Int} # mapped ids 
        nodes::Vector{N}               
        paths::Dict{Int, Dict{Int, Vector{Int}}} 
        edges::Dict{Int, Dict{Int, Int}}         

        function MappedNodeGraph{N}(g::G) where {G <: AbstractGraph, N <: AbstractGraphNode}
            new{N, G}(g, 
                Dict{Int, Int}(), 
                Vector{N}(undef, 0),
                Dict{Int, Dict{Int, Int}}(), 
                Dict{Int, Dict{Int, Int}}(),
            )
        end
    end

    MappedNodeGraph{N}(g=SimpleGraph()) where {N} = MappedNodeGraph{N}(g)
    MappedGraph(::Type{N}) where {N} = MappedNodeGraph{N}(SimpleGraph{Int}())
    MappedDiGraph(::Type{N}) where {N} = MappedNodeGraph{N}(SimpleDiGraph{Int}())

    @inline get_mappedid(g::MappedNodeGraph, node::Int) = g.mid[node]
    @inline get_mappednode(g::MappedNodeGraph, mid::Int) = @inbounds g.nodes[mid]
    @inline get_node(g::MappedNodeGraph, node::Int) = get_mappednode(g, get_mappedid(g, node))

    Base.isempty(g::MappedNodeGraph) = Base.isempty(g.nodes)

    """
        has_vertex(g, node)

    Return true if `node` is contained in the graph `g`.
    """
    @inline has_vertex(g::MappedNodeGraph, node::Int) = haskey(g.mid, node)

    """
        has_path(g, from, to)

    Return true if there is a path between `from` and `to` in the graph `g`.
    """
    function has_path(g::MappedNodeGraph, from::Int, to::Int)
        return has_path(g.graph, get_mappedid(g, from), get_mappedid(g, to))
    end

    """
        add_vertex!(g, node)

    Add `node` to the graph `g`.
    """
    function add_vertex!(g::MappedNodeGraph{T}, node::T) where {T<:AbstractGraphNode}
        nodeid = get_node_id(node)
        has_vertex(g, nodeid) && return 
        
        # insert a new vertex 
        add_vertex!(g.graph)

        # compute mapped id
        mid = nv(g.graph) 
        
        # updates graph
        push!(g.mid, nodeid => mid)
        push!(g.nodes, node)
        nothing
    end

    """
        add_edge!(g::MappedNodeGraph, from::Int, to::Int, [cost])

    Add an edge between `from` and `to` to `g`. 
    Optionally assign a `cost` to the edge.
    """
    function add_edge!(g::MappedNodeGraph{T}, from::Int, to::Int, cost::Int=0) where {T}

        # ensure the two vertexes already exist in the graph 
        if !(has_vertex(g, from) && has_vertex(g, to))
            throw(ErrorException("The vertex $from or $to is not contained in the graph."))
        end

        add_edge!(g.graph, get_mappedid(g, from), get_mappedid(g, to))
        compute_paths!(g)

        edges = get!(g.edges, from, Dict{Int, Int}())
        push!(edges, to => cost)
        nothing 
    end

    """
        compute_paths(g::MappedNodeGraph)

    Compute all possible paths in the graph.
    """
    function compute_paths!(g::MappedNodeGraph{T}) where {T}
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

    """
        get_path(g::MappedNodeGraph, from::Int, to::Int)

    Get the nodes on the path between and including `from` and `to`. Returns an empty array if 
    either `from` or `to` are not a part of `g` or if there is no path between them. 
    """
    function get_path(g::MappedNodeGraph{T}, from::Int, to::Int) where T
        (has_vertex(g, from) && has_vertex(g, to)) || return Int[]
        return g.paths[from][to]
    end

    """
        get_edgecosts(g::MappedNodeGraph, from::Int, to::Int)

    Get all costs assigned to the edges between `from` and `to`. Returns an empty array if 
    either `from` or `to` are not a part of `g` or if there is no path between them.
    """
    function get_edgecosts(g::MappedNodeGraph{T}, from::Int, to::Int) where {T}
        path = get_path(g, from, to)
        isempty(path) && return Int[]
        edges = Vector{Int}(undef, length(path)-1)
        for i in eachindex(edges)
            edges[i] = g.edges[path[i]][path[i+1]]
        end
        return edges 
    end

end