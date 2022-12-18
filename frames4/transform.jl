function Rotation(frame::FrameSystem, from::Union{Int, AbstractFrameAxes}, 
    to::Union{Int, AbstractFrameAxes}, t::Number)

    if from == to 
        #return DCM(1.0I) 
    else 
        fid = from isa Int ? from : axes_id(from)
        tid = to isa Int ? to : axes_id(to)
    
        return _Rotation(frame, t, get_path(frames_axes(frame), fid, tid))
    end

end

@inbounds function _Rotation(frame::FrameSystem, t::Number, path::Vector{Int})

    f1 = get_node(frames_axes(frame), path[1])
    f2 = get_node(frames_axes(frame), path[2]) 

    rot = _Rotation(frame, f1, f2, t)

    for i = 2:length(path)-1
        f1 = f2
        f2 = get_node(frames_axes(frame), path[i+1])
        rot = _Rotation(frame, f1, f2, t) * rot
    end

    return rot 
end

function _Rotation(frames::FrameSystem, from::A, to::A, t::Number) where A
    if from.id == to.parentid
        _Rotation(frames, to, t)
    else
        inv(_Rotation(frames, from, t))
    end
end

function _Rotation(frames::FrameSystem, axes::A, t::Number) where A
    @inbounds if axes.class in (:InertialAxes, :FixedOffsetAxes)
        return axes.R[1]
    else 
        # TODO: implement this
    end
end
