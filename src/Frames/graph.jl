export FrameGraph, FRAMES

import Basic: register!, connect!, find_path

import NodeGraphs: NodeGraph, SimpleGraph, 
        add_edge!, add_vertex!, get_nodes

"""
    FrameGraph

Convenience constructor for `NodeGraph`
"""
function FrameGraph(::Type{N}, ::Type{F}) where {N<:Integer, F}
    NodeGraph{F, N, N}(SimpleGraph{N}())
end
FrameGraph(::Type{F}) where {F} = FrameGraph(Int64, F)


const FRAMES = FrameGraph(AbstractFrame)

"""
    register!(g::NodeGraph, frame::F)

Register a new frame with identifier in the frames graph.
"""
function register!(g::NodeGraph{F, N, G, N}, frame::F) where {N<:Integer, G, F} 
    add_vertex!(g, frame)
end
function register!(frame::F) where {F<:AbstractFrame}
    register!(FRAMES, frame)
end

"""
    connect!(g::NodeGraph{F, N, G, N}, f1::F1, f2::F2) where {N<:Integer, 
        G, F1<:AbstractFrame, F2<:AbstractFrame, F}

Connect two object within the frame graph.
"""
function connect!(g::NodeGraph{F, N, G, N}, f1::F1, f2::F2) where {N<:Integer, 
    G, F1<:AbstractFrame, F2<:AbstractFrame, F}
    add_edge!(g, f1, f2)
end

function connect!(f1::F1, f2::F2) where {F1<:AbstractFrame, F2<:AbstractFrame}
    connect!(FRAMES, f1, f2)
end

"""
    find_path(g::NodeGraph{F, N, G, N}, from::F1, to::F2) where {N<:Integer, 
        G, F1<:AbstractFrame, F2<:AbstractFrame, F}

Find path between two frames in the frame system.
"""
function find_path(g::NodeGraph{F, N, G, N}, from::F1, to::F2) where {N<:Integer, 
    G, F1<:AbstractFrame, F2<:AbstractFrame, F}
    get_nodes(g, from, to)
end

function find_path(f1::F1, f2::F2) where {F1<:AbstractFrame, F2<:AbstractFrame}
    find_path(FRAMES, f1, f2)
end