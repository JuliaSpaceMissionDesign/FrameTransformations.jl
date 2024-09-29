
struct AliasGraph{G,A}
    graph::G
    alias::A
end

"""
    FrameSystem{O, T, S}

A `FrameSystem` instance manages a collection of user-defined `FramePointNode`, 
`FrameAxesNode` and `Direction` objects, enabling computation of arbitrary transformations 
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
    dir::Dict{Symbol,Direction{O,T}}
end

function FrameSystem{O,T,S}() where {O,T,S}
    return FrameSystem{O,T,S}(
        AliasGraph(MappedGraph(FramePointNode{O,T}), Dict{Symbol,Int}()),
        AliasGraph(MappedGraph(FrameAxesNode{O,T}), Dict{Symbol,Int}()),
        Dict()
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

Return the registered points graph.
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