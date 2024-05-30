
"""
    FrameSystem{O, N, S, D}

A `FrameSystem` instance manages a collection of user-defined `FramePointNode` and 
`FrameAxesNode` objects, enabling computation of arbitrary transformations 
between them. It is created by specifying the maximum transformation order `O`, the outputs 
datatype `N` and an `AbstractTimeScale` instance `S`; the parameter `D` is the `length` of 
the output ad shall always be `3O`.

The following transformation orders are accepted: 
- **1**: position 
- **2**: position and velocity 
- **3**: position, velocity and acceleration
- **4**: position, velocity, acceleration and jerk

--- 

    FrameSystem{O, N, S}()

Create a new, empty `FrameSystem` object of order `O`, datatype `N` and timescale `S`.
The parameter `S` can be dropped, in case the default (`BarycentricDynamicalTime`) is used. 
"""
struct FrameSystem{O, N<:Number, S<:AbstractTimeScale, D}
    points::PointsGraph{O, N, D}
    axes::AxesGraph{O, N, D}
end

function FrameSystem{O, N, S}() where {O, N, S}
    D = 3O
    return FrameSystem{O, N, S, D}(
        MappedGraph(FramePointNode{O, N, D}), MappedGraph(FrameAxesNode{O, N, D})
    )
end

@inline FrameSystem{O, N}() where {O, N} = FrameSystem{O, N, BarycentricDynamicalTime}()

function Base.summary(io::IO, ::FrameSystem{O, T, S, D}) where {O,T,S,D}
    return println(io, "FrameSystem{$O, $T, $S, $D}")
end

""" 
    get_order(frames::FrameSystem{O}) where O 

Return the frame system order `O`.
"""
@inline get_order(::FrameSystem{O}) where O = O

""" 
    get_timescale(frames::FrameSystem{O, N, S}) where {O, N, S} 

Return the frame system order timescale `S`.
"""
@inline get_timescale(::FrameSystem{O, N, S}) where {O, N, S} = S 

""" 
    get_points(frames::FrameSystem) 

Return the frame system points graph.
"""
@inline get_points(f::FrameSystem) = f.points

""" 
    get_axes(frames::FrameSystem) 

Return the frame system axes graph.
"""
@inline get_axes(f::FrameSystem) = f.axes

function add_point!(fs::FrameSystem{O, T}, p::FramePointNode{O, T}) where {O,T}
    return add_vertex!(fs.points, p)
end

function add_axes!(fs::FrameSystem{O, T}, ax::FrameAxesNode{O, T}) where {O,T}
    return add_vertex!(fs.axes, ax)
end

""" 
    has_point(frames::FrameSystem, id::Int) 

Check if `id` point is within `frames`.
"""
@inline has_point(f::FrameSystem, id::Int) = has_vertex(get_points(f), id)

""" 
    has_axes(frames::FrameSystem, axesid::Int) 

Check if `axesid` axes is within `frames`.
"""
@inline has_axes(f::FrameSystem, axesid::Int) = has_vertex(get_axes(f), axesid)

# ---
# Formatting & printing 

function _fmt_node(n::FramePointNode)
    return " $(n.name)(id=$(n.id), axesid=$(n.axesid), class=$(n.class))"
end

function _fmt_node(n::FrameAxesNode)
    return " $(n.name)(id=$(n.id), class=$(n.class))"
end

function prettyprint(g::Union{AxesGraph, PointsGraph})
    if !isempty(g.nodes)
        println(_fmt_node(g.nodes[1]))
        _print_frame_graph(g, get_node_id(g.nodes[1]), 2, " ", " │   ")
    end
end

function _print_frame_graph(g, pid::Int, idx::Int, last_prefix::String, prefix::String)
    for i in idx:length(g.nodes)
        if g.nodes[i].parentid == pid
            prefix = (i < length(g.nodes) && g.nodes[i+1].parentid == pid) ? " |" : " └"
            println(last_prefix * prefix * "── " * _fmt_node(g.nodes[i]))
            _print_frame_graph(
                g, get_node_id(g.nodes[i]), i, last_prefix * prefix * "   ", last_prefix * prefix)
        end
    end
end

function Base.show(io::IO, g::FrameSystem{O, N, S, D}) where {O, N, S, D}
    println(
        io, 
        "FrameSystem{$O, $N, $S, $D} with $(length(get_points(g).nodes))" 
        * " points and $(length(get_axes(g).nodes)) axes"
    )
    if !isempty(get_points(g).nodes)
        printstyled(io, "\nPoints: \n"; bold=true)
        prettyprint(get_points(g))
    end
    if !isempty(get_axes(g).nodes)
        printstyled(io, "\nAxes: \n"; bold=true)
        prettyprint(get_axes(g))
    end
end