
struct AliasGraph{G,A}
    graph::G
    alias::A
end

"""
    FrameSystem{O, T, S}

A `FrameSystem` instance manages a collection of user-defined `FramePointNode`, 
`FrameAxesNode` and `DirectionDefinition` objects, enabling computation of arbitrary transformations 
between them. It is created by specifying the maximum transformation order `O`, the outputs 
datatype `T` and an `AbstractTimeScale` instance `S`.

The following transformation orders are accepted: 
- **1**: position 
- **2**: position and velocity 
- **3**: position, velocity and acceleration
- **4**: position, velocity, acceleration and jerk

--- 

    FrameSystem{O, T, S}()

Create a new, empty `FrameSystem` object of order `O`, datatype `T` and timescale `S`.
The parameter `S` can be dropped, in case the default (`BarycentricDynamicalTime`) is used. 
"""
struct FrameSystem{O,T<:Number,S<:AbstractTimeScale}
    points::AliasGraph{PointsGraph{O,T},Dict{Symbol,Int}}
    axes::AliasGraph{AxesGraph{O,T},Dict{Symbol,Int}}
    dir::Dict{Symbol,DirectionDefinition{O,T}}
    _axes_nodes::Dict{Tuple{Int,Int},Vector{FrameAxesNode{O,T}}}
    _points_nodes::Dict{Tuple{Int,Int},Vector{FramePointNode{O,T}}}
end

function FrameSystem{O,T,S}() where {O,T,S}
    return FrameSystem{O,T,S}(
        AliasGraph(MappedGraph(FramePointNode{O,T}), Dict{Symbol,Int}()),
        AliasGraph(MappedGraph(FrameAxesNode{O,T}), Dict{Symbol,Int}()),
        Dict{Symbol,DirectionDefinition{O,T}}(),
        Dict{Tuple{Int,Int},Vector{FrameAxesNode{O,T}}}(),
        Dict{Tuple{Int,Int},Vector{FramePointNode{O,T}}}()
    )
end

@inline FrameSystem{O,T}() where {O,T} = FrameSystem{O,T,BarycentricDynamicalTime}()

function Base.summary(io::IO, ::FrameSystem{O,T,S}) where {O,T,S}
    return println(io, "FrameSystem{$O, $T, $S}")
end

""" 
    order(frames::FrameSystem{O}) where O 

Return the frame system order `O`.
"""
@inline order(::FrameSystem{O}) where {O} = O

""" 
    timescale(frames::FrameSystem{O, T, S}) where {O, T, S} 

Return the frame system order timescale `S`.
"""
@inline timescale(::FrameSystem{O,T,S}) where {O,T,S} = S

""" 
    points_graph(frames::FrameSystem) 

Return the frame system points graph.
"""
@inline points_graph(f::FrameSystem) = f.points.graph

""" 
    axes_graph(frames::FrameSystem) 

Return the frame system axes graph.
"""
@inline axes_graph(f::FrameSystem) = f.axes.graph

"""
    points_alias(f::FrameSystem)

Return the registered points aliases map.
"""
@inline points_alias(f::FrameSystem) = f.points.alias

"""
    axes_alias(f::FrameSystem)

Return the registered axes aliases map.
"""
@inline axes_alias(f::FrameSystem) = f.axes.alias

"""
    directions(f::FrameSystem)

Return the direction dictionary.
"""
@inline directions(f::FrameSystem) = f.dir

"""
    point_id(f::FrameSystem, id)

Get the `id` associate to a point.
"""
@inline point_id(::FrameSystem, id::Int) = id
@inline point_id(f::FrameSystem, name::Symbol) = points_alias(f)[name]

"""
    axes_id(f::FrameSystem, id)

Get the `id` associate to an axes.
"""
@inline axes_id(::FrameSystem, id::Int) = id
@inline axes_id(f::FrameSystem, name::Symbol) = axes_alias(f)[name]

"""
    add_point!(fs::FrameSystem{O, T}, p::FramePointNode{O, T}) where {O,T}

Add point to the frame system.
"""
function add_point!(fs::FrameSystem{O,T}, p::FramePointNode{O,T}) where {O,T}
    push!(fs.points.alias, Pair(p.name, p.id))
    return add_vertex!(fs.points.graph, p)
end

"""
    add_axes!(fs::FrameSystem{O, T}, ax::FrameAxesNode{O, T}) where {O,T}

Add axes to the frame system.
"""
function add_axes!(fs::FrameSystem{O,T}, ax::FrameAxesNode{O,T}) where {O,T}
    push!(fs.axes.alias, Pair(ax.name, ax.id))
    return add_vertex!(fs.axes.graph, ax)
end

""" 
    has_point(frames::FrameSystem, id) 

Check if `id` point is within `frames`.
"""
@inline has_point(f::FrameSystem, id) = has_vertex(points_graph(f), point_id(f, id))

""" 
    has_axes(frames::FrameSystem, ax) 

Check if `ax` axes is within `frames`.
"""
@inline has_axes(f::FrameSystem, ax) = has_vertex(axes_graph(f), axes_id(f, ax))

""" 
    has_axes(frames::FrameSystem, name::Symbol) 

Check if `name` direction is within `frames`.
"""
@inline has_direction(f::FrameSystem, name::Symbol) = haskey(f.dir, name)

"""
    _get_axes_nodes(fr::FrameSystem, fromid::Int, toid::Int)

Return the cached vector of `FrameAxesNode` along the path from `fromid` to `toid`,
or `nothing` if the pair is not in the cache. The cache is lazily built on first access
after any topology change, then reused for all subsequent queries.
"""
function _get_axes_nodes(fr::FrameSystem, fromid::Int, toid::Int)
    if isempty(fr._axes_nodes) && length(axes_graph(fr).nodes) > 1
        _rebuild_axes_cache!(fr)
    end
    return get(fr._axes_nodes, (fromid, toid), nothing)
end

"""
    _get_points_nodes(fr::FrameSystem, fromid::Int, toid::Int)

Return the cached vector of `FramePointNode` along the path from `fromid` to `toid`,
or `nothing` if the pair is not in the cache. The cache is lazily built on first access
after any topology change, then reused for all subsequent queries.
"""
function _get_points_nodes(fr::FrameSystem, fromid::Int, toid::Int)
    if isempty(fr._points_nodes) && length(points_graph(fr).nodes) > 1
        _rebuild_points_cache!(fr)
    end
    return get(fr._points_nodes, (fromid, toid), nothing)
end

# Rebuild the full axes path cache. Cost is O(N²) in the number of registered axes.
# Deferred to first query so that registering N axes costs O(N²) total rather than O(N³).
function _rebuild_axes_cache!(fr::FrameSystem{O,T}) where {O,T}
    empty!(fr._axes_nodes)
    g = axes_graph(fr)
    ids = [n.id for n in g.nodes]
    for fromid in ids, toid in ids
        fromid == toid && continue
        path = get_path(g, fromid, toid)
        isempty(path) && continue
        nodes = Vector{FrameAxesNode{O,T}}(undef, length(path))
        @inbounds for i in eachindex(path)
            nodes[i] = get_mappednode(g, path[i])
        end
        fr._axes_nodes[(fromid, toid)] = nodes
    end
    return nothing
end

# Rebuild the full points path cache. Cost is O(N²) in the number of registered points.
# Deferred to first query so that registering N points costs O(N²) total rather than O(N³).
function _rebuild_points_cache!(fr::FrameSystem{O,T}) where {O,T}
    empty!(fr._points_nodes)
    g = points_graph(fr)
    ids = [n.id for n in g.nodes]
    for fromid in ids, toid in ids
        fromid == toid && continue
        path = get_path(g, fromid, toid)
        isempty(path) && continue
        nodes = Vector{FramePointNode{O,T}}(undef, length(path))
        @inbounds for i in eachindex(path)
            nodes[i] = get_mappednode(g, path[i])
        end
        fr._points_nodes[(fromid, toid)] = nodes
    end
    return nothing
end

# ---
# Formatting & printing 

function _fmt_node(n::FramePointNode)
    return " $(n.name)(id=$(n.id), axesid=$(n.axesid))"
end

function _fmt_node(n::FrameAxesNode)
    return " $(n.name)(id=$(n.id))"
end

function prettyprint(g::Union{AxesGraph,PointsGraph})
    if !isempty(g.nodes)
        println(_fmt_node(g.nodes[1]))
        _print_frame_graph(g, get_node_id(g.nodes[1]), 2, " ", " │   ")
    end
end

function _print_frame_graph(g, pid::Int, idx::Int, last::String, prefix::String)
    for i in idx:length(g.nodes)
        if g.nodes[i].parentid == pid
            prefix = (i < length(g.nodes) && g.nodes[i+1].parentid == pid) ? " |" : " └"
            println(last * prefix * "── " * _fmt_node(g.nodes[i]))
            _print_frame_graph(
                g, get_node_id(g.nodes[i]), i, last * prefix * "   ", last * prefix)
        end
    end
end

function Base.show(io::IO, g::FrameSystem{O,T,S}) where {O,T,S}
    println(
        io,
        "FrameSystem{$O, $T, $S} with $(length(points_graph(g).nodes))"
        *
        " points, $(length(axes_graph(g).nodes)) axes and $(length(g.dir)) directions"
    )
    if !isempty(points_graph(g).nodes)
        printstyled(io, "\nPoints: \n"; bold=true)
        prettyprint(points_graph(g))
    end
    if !isempty(axes_graph(g).nodes)
        printstyled(io, "\nAxes: \n"; bold=true)
        prettyprint(axes_graph(g))
    end
    if !isempty(directions(g))
        printstyled(io, "\nDirections: \n"; bold=true)
        for d in values(directions(g))
            println(" └── $(d.name)(id=$(d.id))")
        end
    end
end