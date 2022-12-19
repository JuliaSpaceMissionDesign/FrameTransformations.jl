# AXES 

function Rotation(frame::FrameSystem, from::Union{Int, AbstractAxes}, 
    to::Union{Int, AbstractAxes}, ep::Number)

    if from == to 
        return DCM(1.0I)
    else 
        fid = from isa Int ? from : frameid(from)
        tid = to isa Int ? to : frameid(to)
    
        return _Rotation(frame, ep, get_path(axes_graph(frame), fid, tid))
    end

end


@inbounds function _Rotation(frame::FrameSystem, ep::Number, path::Vector{Int})

    f1 = get_node(axes_graph(frame), path[1])
    f2 = get_node(axes_graph(frame), path[2]) 

    rot = _rotate(frame, f1, f2, ep)

    for i = 2:length(path)-1
        f1 = f2
        f2 = get_node(axes_graph(frame), path[i+1])

        rot = _rotate(frame, f1, f2, ep)*rot
    end

    return rot 
end

function _rotate(frame::FrameSystem, from::AstroAxes, to::AstroAxes, ep::Number)
    from.id == to.parent ? _rotate(frame, to, ep) : _rotate(frame, from, ep)'
end

function _rotate(frame::FrameSystem, axes::AstroAxes{T}, ep::Number) where T
    @inbounds if axes.class in (0, 1)
        return axes.R[1]
    else
        tid = Threads.threadid()
        if axes.epochs[tid] != ep 
            if axes.class == 2 
                stv = SA{T}[0, 0, 0, 0, 0, 0]
            else
                # uncomment when _vector6 is implemented!
                # stv = _vector6(
                #     get_node_from_id(points_graph(frame), axes.point), ep)

                stv = SA{T}[0, 0, 0, 0, 0, 0]
            end

            axes.R[tid] = axes.fun(ep, stv)
            axes.epochs[tid] = ep 
        end
        return axes.R[tid]
    end
end


# POINTS 

function get_vector3(frame::FrameSystem{T}, from::Union{Int, AbstractPoint}, 
                     to::Union{Int, AbstractPoint}, 
                     axes::Union{Int, AbstractAxes}, 
                     ep::Number) where T
    if from == to 
        return SVector{3, T}(0., 0., 0.)
    else 
        fid = from isa Int ? from : naifid(from)
        tid = to isa Int ? to : naifid(to)
        aid = axes isa Int ? axes : frameid(axes)

        return _get_vector3(frame, ep, aid, 
                    get_path(points_graph(frame), fid, tid))
    end
end 

function _get_vector3(frame::FrameSystem{T}, ep::Number, axes::Int, 
                      path::Vector{Int}) where T

    ps_axes = get_node(points_graph(frame), path[1]).axes
    pe_axes = get_node(points_graph(frame), path[end]).axes

    if axes == ps_axes
        return _get_vector3_backward(frame, ep, path)
    elseif axes == pe_axes
        return _get_vector3_forward(frame, ep, path)
    else 
        # TO BE OPTIMIZED (but probably does not lead to any performance gains)
        return Rotation(frame, pe_axes, axes, ep)*_get_vector3_forward(frame, ep, path)
    end

end

function _get_vector3_forward(frame::FrameSystem, ep::Number, path::Vector{Int})

    p1 = get_node(points_graph(frame), path[1])
    p2 = get_node(points_graph(frame), path[2])

    ax, pos = _vector3(p1, p2, ep)
    for i = 2:length(path)-1
        p1 = p2 
        p2 = get_node(points_graph(frame), path[i+1])

        ax2, pos2 = _vector3(p1, p2, ep)

        # Rotates previous vector to p2's axes
        if ax2 != ax    
            pos = Rotation(frame, ax, ax2, ep)*pos
        end

        # Updates axes and position
        ax = ax2 
        pos += pos2
    end

    return pos 
end

function _get_vector3_backward(frame::FrameSystem, ep::Number, path::Vector{Int})

    p1 = get_node(points_graph(frame), path[end])
    p2 = get_node(points_graph(frame), path[end-1])

    ax, pos = _vector3(p1, p2, ep)
    for i = 2:length(path)-1
        p1 = p2 
        p2 = get_node(points_graph(frame), path[end-i])

        ax2, pos2 = _vector3(p1, p2, ep)

        # Rotates previous vector to p2's axes
        if ax2 != ax 
            pos = Rotation(frame, ax, ax2, ep)*pos
        end 

        # Updates axes and position 
        ax = ax2 
        pos += pos2 

    end 

    return -pos
end

function _vector3(from::AstroPoint, to::AstroPoint, ep::Number)
    if from.NAIFId == to.parent 
        return to.axes, _vector3(to, ep)
    else 
        return from.axes, -_vector3(from, ep)
    end
end

# Missing check on available ephemeris timespan!
# Check is done by CALCEPH which promts to terminal the error, but 
# does not interrupt the flow of the program nor return exit flags
function _vector3(point::AstroPoint, ep::Number)

    @inbounds if point.class == 0 || point.class == 3
        pos = SA[point.stv[1][1], 
                 point.stv[1][2], 
                 point.stv[1][3]]
    else 
        tid = Threads.threadid()
        # if point.epochs[tid] != ep 
            if 0 < point.class < 3 # 1 or 2
                point.epochs[tid] = ep 
                point.fun!(point.stv[tid], ep)
            else 
                # Updatable point has not been updated!
                throw(ErrorException("Updatable Point $(point.name) has not been"*
                    " updatated for epoch $ep."))
            end
        # end 

        pos = SA[point.stv[tid][1], 
                 point.stv[tid][2], 
                 point.stv[tid][3]]
    end

    return pos 
end 
