
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

        # Axes transformations 
        @inline function ($axfun1)(::FrameSystem{<:Any, S1}, from, 
                    to, ::Epoch{S2}) where {S1, S2}
            throw(ArgumentError("Incompatible epoch timescale: expected $S1, found $S2."))
        end

        """
            $($axfun1)(frame::FrameSystem, from, to, ep::Epoch) 

        Compute the rotation that transforms a $(3*$order)-elements state vector from one 
        specified set of axes to another at a given epoch. 

        ### Inputs 
        - `frame` -- The `FrameSystem` container object 
        - `from` -- ID or instance of the axes to transform from 
        - `to` -- ID or instance of the axes to transform to 
        - `ep` -- `Epoch` of the rotation. Its timescale must match that of the frame system. 

        ### Output
        A `Rotation` object of order $($order).
        """
        @inline function ($axfun1)(frame::FrameSystem{<:Any, S}, from, to, 
                            ep::Epoch{S}) where S
            $(axfun1)(frame, from, to, j2000(ep))
        end

        """
            $($axfun1)(frame::FrameSystem, from, to, t::Number)
        
        Compute the rotation that transforms a $(3*$order)-elements state vector from one 
        specified set of axes to another at a given time, expressed in days since [`J2000`](@ref) 
        if ephemerides are used. 
        """
        function ($axfun1)(frame::FrameSystem{T}, from, to, t::Number) where T
            from == to && return Rotation{$order}(T(1)I)
            $(axfun2)(frame, t, get_path(frames_axes(frame), 
                                          axes_alias(from), axes_alias(to)))
        end

        # Low-level function to parse a path of axes and chain their rotations 
        @inbounds function ($axfun2)(frame::FrameSystem, t::Number, path::Vector{Int})
            f1 = get_mappednode(frames_axes(frame), path[1])
            f2 = get_mappednode(frames_axes(frame), path[2])

            rot = $(axfun2)(frame, f1, f2, t)

            for i = 2:length(path)-1
                f1 = f2
                f2 = get_mappednode(frames_axes(frame), path[i+1])

                rot = $(axfun2)(frame, f1, f2, t)*rot
            end

            return rot 
        end

        # Low-level function to compute the rotation between two axes
        @inline function ($axfun2)(frame::FrameSystem, from::FrameAxesNode, 
                            to::FrameAxesNode, t::Number)
            from.id == to.parentid ? $(axfun2)(frame, to, t) : inv($(axfun2)(frame, from, t))
        end

        # Low-level function to compute the rotation matrix of a specific set of axes 
        function ($axfun2)(frame::FrameSystem, axes::FrameAxesNode{T}, t::Number) where T 
            @inbounds if axes.class in (:InertialAxes, :FixedOffsetAxes)
                return $order < 3 ? Rotation{$order}(axes.R[1]) : axes.R[1]
            else 
                tid = Threads.threadid()
                if axes.epochs[tid] != t || axes.nzo[tid] < $order
                    if axes.class == :RotatingAxes 
                        stv = @SVector zeros(T, 3*$order)
                        axes.R[tid] = axes.$(f)(t, stv, stv)
                    else 
                        axes.R[tid] = axes.$(f)(t, 
                            $(compfun)(frame, axes.comp.v1, axes.parentid, t), 
                            $(compfun)(frame, axes.comp.v2, axes.parentid, t))
                    end 
                    axes.epochs[tid] = t 
                    axes.nzo[tid] = $order

                end

                return $order < 3 ? Rotation{$order}(axes.R[tid]) : axes.R[tid]
            end
        end


        # Point transformations
        
        """
            $($pfun1)(frame::FrameSystem, from, to, axes, ep::Epoch) 

        Compute $(3*$order)-elements state vector of a target point relative to 
        an observing point, in a given set of axes, at the desired epoch.

        ### Inputs 
        - `frame` -- The `FrameSystem` container object 
        - `from` -- ID or instance of the observing point
        - `to` -- ID or instance of the target point 
        - `axes` -- ID or instance of the output state vector axes 
        - `ep` -- `Epoch` of the observer. Its timescale must match that of the frame system. 

        """
        @inline function ($pfun1)(::FrameSystem{<:Any, S1}, from, 
                            to, axes, ::Epoch{S2}) where {S1, S2}
            throw(ArgumentError("Incompatible epoch timescale: expected $S1, found $S2."))
        end

        @inline function ($pfun1)(frame::FrameSystem{<:Any, S}, from, 
                            to, axes, ep::Epoch{S}) where S
            $(pfun1)(frame, from, to, axes, j2000(ep))
        end

        """
            $($pfun1)(frame::FrameSystem, from, to, axes, t::Number) 

        Compute $(3*$order)-elements state vector of a target point relative to 
        an observing point, in a given set of axes, at the desired time expressed in 
        days since [`J2000`](@ref) if ephemerides are used. 
        """
        function ($pfun1)(frame::FrameSystem{T}, from, to, axes, t::Number) where T

            from == to && return @SVector zeros(T, 3*$order)
            $(pfun2)(frame, t, axes_alias(axes), get_path(frames_points(frame), 
                     point_alias(from), point_alias(to)))

        end

        # Low-level function to select the most-convenient parsing direction of a given 
        # path of points depending on the input axes
        function ($pfun2)(frame::FrameSystem, t::Number, axesid::Int, path::Vector{Int})

            @inbounds ps_axes = get_mappednode(frames_points(frame), path[1]).axesid
            @inbounds pe_axes = get_mappednode(frames_points(frame), path[end]).axesid

            if axesid == ps_axes
                return $(vbwd)(frame, t, path)
            elseif axesid == pe_axes
                return $(vfwd)(frame, t, path)
            else 
                # TO BE OPTIMIZED (but probably does not lead to any performance gains)
                return $(axfun1)(frame, pe_axes, axesid, t)*$(vfwd)(frame, t, path)
            end
        end

        # Low-level function to chain point translations in a forward direction 
        @inbounds function ($vfwd)(frame::FrameSystem, t::Number, path::Vector{Int}) 
            p1 = get_mappednode(frames_points(frame), path[1])
            p2 = get_mappednode(frames_points(frame), path[2])
        
            axid, stv = ($pfun2)(p1, p2, t)
            for i = 2:length(path)-1
                p1 = p2 
                p2 = get_mappednode(frames_points(frame), path[i+1])
        
                ax2id, stv2 = ($pfun2)(p1, p2, t)
        
                # Rotates previous vector to p2's axes
                if ax2id != axid    
                    stv = ($axfun1)(frame, axid, ax2id, t)*stv
                end
        
                # Updates axes and position
                axid = ax2id 
                stv += stv2
            end
        
            return stv 
        end

        # Low-level function to chain point translations in a backward direction 
        @inbounds function ($vbwd)(frame::FrameSystem, t::Number, path::Vector{Int}) 
            p1 = get_mappednode(frames_points(frame), path[end])
            p2 = get_mappednode(frames_points(frame), path[end-1])
        
            axid, stv = ($pfun2)(p1, p2, t)
            for i = 2:length(path)-1
                p1 = p2 
                p2 = get_mappednode(frames_points(frame), path[end-i])
        
                ax2id, stv2 = ($pfun2)(p1, p2, t)
        
                # Rotates previous vector to p2's axes
                if ax2id != axid    
                    stv = ($axfun1)(frame, axid, ax2id, t)*stv
                end
        
                # Updates axes and position
                axid = ax2id 
                stv += stv2
            end
        
            return -stv 
        end

        # Low-level function to compute the translation between two points
        function ($pfun2)(from::FramePointNode, to::FramePointNode, t::Number)
            if from.NAIFId == to.parentid 
                return to.axesid, $(pfun2)(to, t)
            else 
                return from.axesid, -$(pfun2)(from, t)
            end
        end

        # Low-level function to compute the translation of a given point
        @inbounds function ($pfun2)(point::FramePointNode{T}, t::Number) where T 
            if point.class in (:RootPoint, :FixedPoint)
                return SA[point.stv[1].data[1:3*$order]...]
            else
                tid = Threads.threadid()
                if point.epochs[tid] != t || point.nzo[tid] < $order 
                    if point.class == :UpdatablePoint
                        # Updatable point has not been updated! 
                        throw(ErrorException(
                            "UpdatablePoint with NAIFId = $(point.NAIFId) has not been "*
                            "updated at time $t for order $($order)"))
                    else 
                        point.$(f)(point.stv[tid], t)
                    end

                    point.epochs[tid] = t 
                    point.nzo[tid] = $order 
                end

                return SA[point.stv[tid].data[1:3*$order]...]
            end
        end

    end
end


""" 
    update_point!(frames, point, stv::AbstractVector, epoch::Epoch)
"""
function update_point!(frames::FrameSystem{T, S}, point, stv::AbstractVector{T}, 
            epoch::Epoch{S}) where {T, S}
    update_point!(frames, point, stv, j2000(epoch))
end

""" 
    update_point!(frames::FrameSystem{T}, point, stv::AbstractVector{T}, time::T) where T

Update the state vector of `point` at the input `time` in `frames`. The only 
accepted length for the input vector `stv` are 3, 6 or 9. The order is automatically inferred 
from the vector length.

### Examples 
```jldoctest
julia> FRAMES = FrameSystem{Float64}();
  
julia> @axes ICRF 1  

julia> add_axes_inertial!(FRAMES, ICRF)

julia> @point Origin 0

julia> @point Satellite 1 

julia> add_point_root!(FRAMES, Origin, ICRF)

julia> add_point_updatable!(FRAMES, Satellite, Origin, ICRF)

julia> y = [10000., 200., 300.]

julia> update_point!(FRAMES, Satellite, y, 0.1)

julia> get_vector3(FRAMES, Origin, Satellite, ICRF, 0.1)
3-element SVector{3, Float64} with indices SOneTo(3):
 10000.0
   200.0
   300.0

julia> get_vector3(FRAMES, Origin, Satellite, ICRF, 0.2)
ERROR: UpdatablePoint with NAIFId = 1 has not been updated at time 0.2 for order 1

julia> get_vector6(FRAMES, Origin, Satellite, ICRF, 0.1)
ERROR: UpdatablePoint with NAIFId = 1 has not been updated at time 0.2 for order 2
```
### See also 
See also [`add_point_updatable!`](@ref)
"""
function update_point!(frames::FrameSystem{T}, point, stv::AbstractVector{T}, 
            time::T) where T

    NAIFId = point_alias(point) 

    !has_point(frames, NAIFId) && throw(ErrorException(
        "Point with NAIF ID $NAIFId is not contained in the frame system."))

    ne = length(stv) 

    if !(ne in (3, 6, 9)) 
        throw(ArgumentError(
            "wrong state vector length: expected either 3, 6 or 9, found $ne"))
    end 

    pnt = get_node(frames_points(frames), point_alias(point))

    id = Threads.threadid()
    @inbounds pnt.epochs[id] = time
    @inbounds pnt.nzo[id] = ne ÷ 3 
    @inbounds @views pnt.stv[id][1:ne] .= stv[1:ne]
    nothing 

end


""" 
    _get_comp_axes_vector3(frame, v, axesid, t)

Compute a 3-elements vector in the desired axes at the given time 
between two points of the frame system 

The returned vector depends on the order in `v` as follows: 

- **1**: position
- **2**: velocity
- **3**: acceleration 

"""
function _get_comp_axes_vector3(frame::FrameSystem, v::ComputableAxesVector, 
                        axesid::Int, t::Number)
    
    @inbounds if v.order == 1 
        return get_vector3(frame, v.from, v.to, axesid, t)        
    elseif v.order == 2 
        stv = get_vector6(frame, v.from, v.to, axesid, t)
        return SA[stv[4], stv[5], stv[6]]
    else 
        stv = get_vector9(frame, v.from, v.to, axesid, t)
        return SA[stv[7], stv[8], stv[9]]
    end

end

""" 
    _get_comp_axes_vector6(frame, v, axesid, t)

Compute a 6-elements vector in the desired axes at the given time 
between two points of the frame system 

The returned vector depends on the order in `v` as follows: 

- **1**: position, velocity
- **2**: velocity, acceleration 
- **3**: acceleration, jerk

### Notes 
This function only returns vectors up to order 2, because the 
frame system currently is uncapable of computing the jerk.
"""
 function _get_comp_axes_vector6(frame::FrameSystem, v::ComputableAxesVector, 
                        axesid::Int, t::Number)
    
    @inbounds if v.order == 1 
        return get_vector6(frame, v.from, v.to, axesid, t)        
    elseif v.order == 2 
        stv = get_vector9(frame, v.from, v.to, axesid, t)
        return SA[stv[4], stv[5], stv[6], stv[7], stv[8], stv[9]]
    end

    throw(ErrorException(
        "Unable to compute a vector of order 4 (jerk). The maximum available order is 3."))
end

""" 
    _get_comp_axes_vector9(frame, v, axesid, t)

Compute a 9-elements vector in the desired axes at the given time 
between two points of the frame system 

The returned vector depends on the order in `v` as follows: 

- **1**: position, velocity, acceleration
- **2**: velocity, acceleration, jerk
- **3**: acceleration, jerk, jounce

### Notes 
This function only works for vectors of order 1, because the 
frame system currently is uncapable of computing the jerk and jounce.
"""
function _get_comp_axes_vector9(frame::FrameSystem, v::ComputableAxesVector, 
            axesid::Int, t::Number)
    
    if v.order == 1 
        return get_vector9(frame, v.from, v.to, axesid, t)        
    end

    throw(ErrorException("Unable to compute a vector of order $(2+v.order). "*
        "The maximum available order is 3."))
end


