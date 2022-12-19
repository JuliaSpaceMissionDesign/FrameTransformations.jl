# AXES 

for (order, f, axfun1, axfun2, pfun1, pfun2, compfun, vfwd, vbwd) in zip(
        (1, 2, 3), (:f, :δf, :δ²f),  
        (:get_rotation3, :get_rotation6, :get_rotation9), 
        (:_compute_rot3, :_compute_rot6, :_compute_rot9),
        (:get_vector3, :get_vector6, :get_vector9), 
        (:_compute_vector3, :_compute_vector6, :_compute_vector9),
        (:_get_comp_axes_vector3, :_get_comp_axes_vector6, :_get_comp_axes_vector9),
        (:_get_vector3_forward, :_get_vector6_forward, :_get_vector9_forward),
        (:_get_vector3_backwards, :_get_vector6_backwards, :_get_vector9_backwards)
    )

    @eval begin 
        function ($axfun1)(frame::FrameSystem{T, <:Any}, from, to, ep::Number) where T
            from == to && return Rotation{$order}(T(1)I)
            $(axfun2)(frame, ep, get_path(frames_axes(frame), axes_alias(from), axes_alias(to)))
        end

        function ($axfun2)(frame::FrameSystem, ep::Number, path::Vector{Int})
            @inbounds f1 = get_mappednode(frames_axes(frame), path[1])
            @inbounds f2 = get_mappednode(frames_axes(frame), path[2])

            rot = $(axfun2)(frame, f1, f2, ep)

            @inbounds for i = 2:length(path)-1
                f1 = f2
                f2 = get_mappednode(frames_axes(frame), path[i+1])

                rot = $(axfun2)(frame, f1, f2, ep)*rot
            end

            return rot 
        end

        @inline function ($axfun2)(frame::FrameSystem, from::FrameAxesNode, 
                            to::FrameAxesNode, ep::Number)
            from.id == to.parentid ? $(axfun2)(frame, to, ep) : inv($(axfun2)(frame, from, ep))
        end

        function ($axfun2)(frame::FrameSystem, axes::FrameAxesNode{T}, ep::Number) where T 
            @inbounds if axes.class in (:InertialAxes, :FixedOffsetAxes)
                return $order < 3 ? Rotation{$order}(axes.R[1]) : axes.R[1]
            else 
                tid = Threads.threadid()
                if axes.epochs[tid] != ep || axes.nzo[tid] < $order
                    if axes.class == :RotatingAxes 
                        stv = @SVector zeros(T, 3*$order)
                        axes.R[tid] = axes.$(f)(ep, stv, stv)
                    else 
                        axes.R[tid] = axes.$(f)(ep, 
                            $(compfun)(frame, axes.comp.v1, axes.parentid, ep), 
                            $(compfun)(frame, axes.comp.v2, axes.parentid, ep))
                    end 
                    axes.epochs[tid] = ep 
                    axes.nzo[tid] = $order

                end

                return $order < 3 ? Rotation{$order}(axes.R[tid]) : axes.R[tid]
            end
        end

        # Point 

        function ($pfun1)(frame::FrameSystem{T, <:Any}, from, to, axes, ep::Number) where T
            from == to && return @SVector zeros(T, 3*$order)
            $(pfun2)(frame, ep, axes_alias(axes), get_path(frames_points(frame), point_alias(from), point_alias(to)))
        end


        function ($pfun2)(frame::FrameSystem, ep::Number, axesid::Int, path::Vector{Int})

            @inbounds ps_axes = get_mappednode(frames_points(frame), path[1]).axesid
            @inbounds pe_axes = get_mappednode(frames_points(frame), path[end]).axesid

            if axesid == ps_axes
                return $(vbwd)(frame, ep, path)
            elseif axesid == pe_axes
                return $(vfwd)(frame, ep, path)
            else 
                # TO BE OPTIMIZED (but probably does not lead to any performance gains)
                return $(axfun1)(frame, pe_axes, axesid, ep)*$(vfwd)(frame, ep, path)
            end
        end

        function ($vfwd)(frame::FrameSystem, ep::Number, path::Vector{Int}) 
            p1 = get_mappednode(frames_points(frame), path[1])
            p2 = get_mappednode(frames_points(frame), path[2])
        
            axid, stv = ($pfun2)(p1, p2, ep)
            for i = 2:length(path)-1
                p1 = p2 
                p2 = get_mappednode(frames_points(frame), path[i+1])
        
                ax2id, stv2 = ($pfun2)(p1, p2, ep)
        
                # Rotates previous vector to p2's axes
                if ax2id != axid    
                    stv = ($axfun1)(frame, axid, ax2id, ep)*stv
                end
        
                # Updates axes and position
                axid = ax2id 
                stv += stv2
            end
        
            return stv 
        end

        function ($vbwd)(frame::FrameSystem, ep::Number, path::Vector{Int}) 
            p1 = get_mappednode(frames_points(frame), path[end])
            p2 = get_mappednode(frames_points(frame), path[end-1])
        
            axid, stv = ($pfun2)(p1, p2, ep)
            for i = 2:length(path)-1
                p1 = p2 
                p2 = get_mappednode(frames_points(frame), path[end-i])
        
                ax2id, stv2 = ($pfun2)(p1, p2, ep)
        
                # Rotates previous vector to p2's axes
                if ax2id != axid    
                    stv = ($axfun1)(frame, axid, ax2id, ep)*stv
                end
        
                # Updates axes and position
                axid = ax2id 
                stv += stv2
            end
        
            return -stv 
        end

        function ($pfun2)(from::FramePointNode, to::FramePointNode, ep::Number)
            if from.NAIFId == to.parentid 
                return to.axesid, $(pfun2)(to, ep)
            else 
                return from.axesid, -$(pfun2)(from, ep)
            end
        end

        function ($pfun2)(point::FramePointNode{T}, ep::Number) where T 
            @inbounds if point.class in (:RootPoint, :FixedPoint)
                return SA[point.stv[1].data[1:3*$order]...]
            else
                tid = Threads.threadid()
                # if point.epochs[tid] != ep || point.nzo[tid] < $order 
                    if point.class == :UpdatablePoint
                        # Updatable point has not been updated! 
                        throw(ErrorException(
                            "UpdatablePoint with NAIFId = $(point.NAIFId) has not been "*
                            "updated for epoch $ep at order $($order)"))
                    else 
                        point.$(f)(point.stv[tid], ep)
                    end

                    point.epochs[tid] = ep 
                    point.nzo[tid] = $order 
                # end

                return SA[point.stv[tid].data[1:3*$order]...]
            end
        end

    end
end


function _get_comp_axes_vector3(frame::FrameSystem, v::ComputableAxesVector, axesid::Int, ep::Number)
    
    if v.order == 1 
        return get_vector_3(frame, v.to, v.from, axesid, ep)        
    elseif v.order == 2 
        stv = get_vector_6(frame, v.to, v.from, axesid, ep)
        return SA[stv[4], stv[5], stv[6]]
    else 
        stv = get_vector_9(frame, v.to, v.from, axesid, ep)
        return SA[stv[7], stv[8], stv[9]]
    end

end


function _get_comp_axes_vector6(frame::FrameSystem, v::ComputableAxesVector, axesid::Int, ep::Number)
    
    if v.order == 1 
        return get_vector_6(frame, v.to, v.from, axesid, ep)        
    elseif v.order == 2 
        stv = get_vector_9(frame, v.to, v.from, axesid, ep)
        return SA[stv[4], stv[5], stv[6], stv[7], stv[8], stv[9]]
    end

    throw(ErrorException(
        "Unable to compute a vector of order 4 (jerk). The maximum available order is 3."))
end

function _get_comp_axes_vector9(frame::FrameSystem, v::ComputableAxesVector, axesid::Int, ep::Number)
    
    if v.order == 1 
        return get_vector_9(frame, v.to, v.from, axesid, ep)        
    end

    throw(ErrorException(
        "Unable to compute a vector of order $(2+v.order). The maximum available order is 3."))
end


