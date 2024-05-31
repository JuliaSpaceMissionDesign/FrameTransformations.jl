
"""
    FrameSystem{O, N, S, D}

A `FrameSystem` instance manages a collection of user-defined `FramePointNode`, 
`FrameAxesNode` and `Direction` objects, enabling computation of arbitrary transformations 
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
    dir::Dict{Symbol, Direction{O, N, D}}

    points_map::Dict{Symbol, Int}
    axes_map::Dict{Symbol, Int}
end

function FrameSystem{O, N, S}() where {O, N, S}
    D = 3O
    return FrameSystem{O, N, S, D}(
        MappedGraph(FramePointNode{O, N, D}), MappedGraph(FrameAxesNode{O, N, D}), Dict(),
        Dict(), Dict()
    )
end

@inline FrameSystem{O, N}() where {O, N} = FrameSystem{O, N, BarycentricDynamicalTime}()

function Base.summary(io::IO, ::FrameSystem{O, T, S, D}) where {O,T,S,D}
    return println(io, "FrameSystem{$O, $T, $S, $D}")
end

""" 
    order(frames::FrameSystem{O}) where O 

Return the frame system order `O`.
"""
@inline order(::FrameSystem{O}) where O = O

""" 
    timescale(frames::FrameSystem{O, N, S}) where {O, N, S} 

Return the frame system order timescale `S`.
"""
@inline timescale(::FrameSystem{O, N, S}) where {O, N, S} = S 

""" 
    points_graph(frames::FrameSystem) 

Return the frame system points graph.
"""
@inline points_graph(f::FrameSystem) = f.points

""" 
    axes_graph(frames::FrameSystem) 

Return the frame system axes graph.
"""
@inline axes_graph(f::FrameSystem) = f.axes

"""
    directions_map(f::FrameSystem)

Return the direction dictionary.
"""
@inline directions_map(f::FrameSystem) = f.dir

"""
    axes(f::FrameSystem)

Return the registered axes names/ids map.
"""
@inline Base.axes(f::FrameSystem) = f.axes_map

"""
    points(f::FrameSystem)

Return the registered points names/ids map.
"""
@inline points(f::FrameSystem) = f.points_map

"""
    directions(f::FrameSystem)

Return the registered directions names.
"""
@inline directions(f::FrameSystem) = keys(directions_map(f))

"""
    point_id(f::FrameSystem, id)

Get the `id` associate to a point.
"""
@inline point_id(::FrameSystem, id::Int) = id 
@inline point_id(f::FrameSystem, name::Symbol) = points(f)[name]

"""
    axes_id(f::FrameSystem, id)

Get the `id` associate to an axes.
"""
@inline axes_id(::FrameSystem, id::Int) = id
@inline axes_id(f::FrameSystem, name::Symbol) = axes(f)[name]

function add_point!(fs::FrameSystem{O, T}, p::FramePointNode{O, T}) where {O,T}
    push!(fs.points_map, Pair(p.name, p.id))
    return add_vertex!(fs.points, p)
end

function add_axes!(fs::FrameSystem{O, T}, ax::FrameAxesNode{O, T}) where {O,T}
    push!(fs.axes_map, Pair(ax.name, ax.id))
    return add_vertex!(fs.axes, ax)
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
        "FrameSystem{$O, $N, $S, $D} with $(length(points_graph(g).nodes))" 
        * " points, $(length(axes_graph(g).nodes)) axes and $(length(g.dir)) directions"
    )
    if !isempty(points_graph(g).nodes)
        printstyled(io, "\nPoints: \n"; bold=true)
        prettyprint(points_graph(g))
    end
    if !isempty(axes_graph(g).nodes)
        printstyled(io, "\nAxes: \n"; bold=true)
        prettyprint(axes_graph(g))
    end
    if !isempty(directions_map(g))
        printstyled(io, "\nDirections: \n"; bold=true)
        for d in values(directions_map(g))
            println(" └── $(d.name)(id=$(d.id))")
        end
    end
end