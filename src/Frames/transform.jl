export vector3, vector6, vector9, vector12,
       rotation3, rotation6, rotation9, rotation12,
       update_point!

for (order, axfun1, axfun2, pfun1, pfun2, compfun, vfwd, vbwd) in zip(
        (1, 2, 3, 4), 

        (:rotation3, :rotation6, :rotation9, :rotation12), 
        (:_compute_rot3, :_compute_rot6, :_compute_rot9, :_compute_rot12),
        
        (:vector3, :vector6, :vector9, :vector12), 
        (:_compute_vector3, :_compute_vector6, :_compute_vector9, :_compute_vector12),

        (:_get_comp_axes_vector3, :_get_comp_axes_vector6, 
         :_get_comp_axes_vector9, :_get_comp_axes_vector12),

        (:_vector3_forward, :_vector6_forward, 
         :_vector9_forward, :_vector12_forward),

        (:_vector3_backwards, :_vector6_backwards, 
         :_vector9_backwards, :_vector12_backwards)
    )

    @eval begin 

        # ----------------------------------
        # AXES TRANSFORMATIONS 
        # ----------------------------------
        @inline function ($axfun1)(::FrameSystem{<:Any, <:Any, S1}, from, 
                    to, ::Epoch{S2}) where {S1, S2}
            throw(ArgumentError("Incompatible epoch timescale: expected $S1, found $S2."))
        end

        """
            $($axfun1)(frame::FrameSystem, from, to, ep::Epoch) 

        Compute the rotation that transforms a $(3*$order)-elements state vector from one 
        specified set of axes to another at a given epoch. 

        Requires a frame system of order ≥ $($order).

        ### Inputs 
        - `frame` -- The `FrameSystem` container object 
        - `from` -- ID or instance of the axes to transform from 
        - `to` -- ID or instance of the axes to transform to 
        - `ep` -- `Epoch` of the rotation. Its timescale must match that of the frame system. 

        ### Output
        A `Rotation` object of order $($order).
        """
        @inline function ($axfun1)(frame::FrameSystem{<:Any, <:Any, S}, from, to, 
                            ep::Epoch{S}) where S
            $(axfun1)(frame, from, to, Tempo.j2000s(ep))
        end

        """
            $($axfun1)(frame::FrameSystem, from, to, t::Number)
        
        Compute the rotation that transforms a $(3*$order)-elements state vector from one 
        specified set of axes to another at a given time `t`, expressed in seconds since 
        [`J2000`](@ref) if ephemerides are used. 
        """
        function ($axfun1)(frame::FrameSystem{O, T}, from, to, t::Number) where {O, T}
            
            if O < $order 
                throw(ErrorException("Insufficient frame system order: "*
                    "transformation requires at least order $($order)."))
            end

            # Retrieve aliased axes ids
            idfrom = axes_alias(from)
            idto = axes_alias(to) 

            idfrom == idto && return Rotation{$order}(T(1)I)
            $(axfun2)(frame, t, get_path(frames_axes(frame), idfrom, idto))
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
        function ($axfun2)(frame::FrameSystem{O, T}, axes::FrameAxesNode{O, T}, 
                    t::Number) where {O, T} 
                    
            @inbounds if axes.class in (:InertialAxes, :FixedOffsetAxes)
                return $order < O ? Rotation{$order}(axes.R[1]) : axes.R[1]
            else 
                tid = Threads.threadid()
                if axes.epochs[tid] != t || axes.nzo[tid] < $order

                    if axes.class in (:RotatingAxes, :ProjectedAxes)
                        stv = @SVector zeros(T, 3O)
                        axes.R[tid] = axes.f[$order](t, stv, stv)

                    elseif axes.class == :EphemerisAxes 
                        stv = @SVector zeros(T, 3O)

                        # Retrieves libration angles for ephemeris kernels 
                        ephem_orient_order!(axes.angles[tid], frame.eph, DJ2000, 
                                            t/DAY2SEC, axes.id, $order-1)
                        
                        # Compute rotation matrix
                        axes.R[tid] = axes.f[$order](t, SA[axes.angles[tid]...], stv)

                    else 
                        axes.R[tid] = axes.f[$order](t, 
                            $(compfun)(frame, axes.comp.v1, axes.parentid, t), 
                            $(compfun)(frame, axes.comp.v2, axes.parentid, t))
                    end 

                    axes.epochs[tid] = t 
                    axes.nzo[tid] = $order
                end

                return $order < O ? Rotation{$order}(axes.R[tid]) : axes.R[tid]
            end
        end


        # ----------------------------------
        # POINTS TRANSFORMATIONS 
        # ----------------------------------
        """
            $($pfun1)(frame::FrameSystem, from, to, axes, ep::Epoch) 

        Compute $(3*$order)-elements state vector of a target point relative to 
        an observing point, in a given set of axes, at the desired epoch `ep`.

        Requires a frame system of order ≥ $($order).

        ### Inputs 
        - `frame` -- The `FrameSystem` container object 
        - `from` -- ID or instance of the observing point
        - `to` -- ID or instance of the target point 
        - `axes` -- ID or instance of the output state vector axes 
        - `ep` -- `Epoch` of the observer. Its timescale must match that of the frame system. 

        """
        @inline function ($pfun1)(::FrameSystem{<:Any, <:Any, S1}, from, 
            to, axes, ::Epoch{S2}) where {S1, S2}
            throw(ArgumentError("Incompatible epoch timescale: expected $S1, found $S2."))
        end

        @inline function ($pfun1)(frame::FrameSystem{<:Any, <:Any, S}, from, 
                    to, axes, ep::Epoch{S}) where S
            $(pfun1)(frame, from, to, axes, Tempo.j2000s(ep))
        end

        """
            $($pfun1)(frame, from, to, axes, t::Number) 

        Compute $(3*$order)-elements state vector of a target point relative to 
        an observing point, in a given set of axes, at the desired time `t` expressed in 
        seconds since [`J2000`](@ref) if ephemerides are used. 

        """
        function ($pfun1)(frame::FrameSystem{O, T}, from, to, axes, t::Number) where {O, T}

            if O < $order 
                throw(ErrorException("Insufficient frame system order: "*
                    "transformation requires at least order $($order)."))
            end

            # Retrieve aliased point ids
            idfrom = point_alias(from)
            idto = point_alias(to) 

            from == to && return @SVector zeros(T, 3*$order)
            $(pfun2)(frame, t, axes_alias(axes), get_path(frames_points(frame), 
                        idfrom, idto))

        end

        # Low-level function to select the most-convenient parsing direction of a given 
        # path of points depending on the input axes
        function ($pfun2)(frame::FrameSystem, t::Number, axesid::Int, path::Vector{Int})

            @inbounds p1 = get_mappednode(frames_points(frame), path[1])
            @inbounds p2 = get_mappednode(frames_points(frame), path[end])
            
            if length(path) == 2
                # This handles all the cases where you don't need to chain any transformations
                ax2id, stv = ($pfun2)(p1, p2, t)
                if ax2id == axesid 
                    return stv
                else 
                    return $(axfun1)(frame, ax2id, axesid, t)*stv 
                end
            
            elseif axesid == p1.axesid
                return $(vbwd)(frame, t, path, p2)
            
            elseif axesid == p2.axesid
                return $(vfwd)(frame, t, path, p1)
            
            else 
                # Optimising this transformation would probably demand a significant 
                # portion of time with respect to the time required by the whole transformation
                return $(axfun1)(frame, p2.axesid, axesid, t)*$(vfwd)(frame, t, path, p1)
            end
        end

        # Low-level function to chain point translations in a forward direction 
        @inbounds function ($vfwd)(frame::FrameSystem, t::Number, path::Vector{Int}, 
                    p1::FramePointNode)

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
        @inbounds function ($vbwd)(frame::FrameSystem, t::Number, path::Vector{Int}, 
                    p1::FramePointNode) 

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
        @inbounds function ($pfun2)(point::FramePointNode{O, T}, t::Number) where {O, T}  
            if point.class in (:RootPoint, :FixedPoint)
                return SA[point.stv[1].data[1:3*$order]...]
            else
                tid = Threads.threadid()

                if point.epochs[tid] != t || point.nzo[tid] < $order 
                    if point.class == :UpdatablePoint
                        # Updatable point has not been updated! 
                        throw(ErrorException(
                            "UpdatablePoint with NAIFId $(point.NAIFId) has not been "*
                            "updated at time $t for order $($order)"))
                    else 
                        point.f[$order](point.stv[tid], t)
                    end

                    point.epochs[tid] = t 
                    point.nzo[tid] = $order 
                end

                return SA[point.stv[tid].data[1:3*$order]...]
            end
        end

    end
end


# ---------------------------------------------
# LIGHT TIME and STELLAR ABERRATION CORRECTIONS 
# ---------------------------------------------

for (order, pfun, ltcorr) in zip(
        (1, 2), (:vector3, :vector6), (:light_time_corr3, :light_time_corr6)
    )

    @eval begin 

        @inline function ($pfun)(::FrameSystem{<:Any, <:Any, S1}, from, to, axes, 
            ::Epoch{S2}, ::AbstractLightTimeCorrection, dir::Int; kwargs...) where {S1, S2}
            throw(ArgumentError("Incompatible epoch timescale: expected $S1, found $S2."))
        end


        """
            $($pfun)(frame, from, to, axes, ep::Epoch, ltcorr, dir; <keyword arguments>) 

        Compute a light-time corrected $(3*$order)-elements state vector of a target point 
        relative to an observing point, in a given set of axes, at the desired epoch `ep`, 
        using the aberration flag `ltcorr`, which may be any of the following `AbstractLightTimeCorrection`:
        
        - **LightTimeCorrection**: it applies the one-way light time (planetary aberration) 
            correction, using a Newtonian formulation. 

        - **StellarAberrationCorrection**: it applies the one-way light time and stellar 
            aberration corrections using a Newtonian fromulation. It modifies the vector 
            obtained with the `LightTimeCorrection` option to account for the observer velocity 
            with respect to the Solar System Barycenter. 
        
        The integer argument `dir` is used to specify correction direction, as follows:
            
        - **-1**: for **Reception**, in which photons depart from the target's location at the 
            light-time corrected epoch `ep-lt` and arrive at the observer's location at `ep`.

        - **+1**: for **Transmission**, in which photons depart from the observer's location at
            `ep` and arrive at the target's location at the light-time corrected epoch `ep+lt`.
    
        ### Keyword Arguments 
                
        - `iters::Int=1`: the number of iterations used to find the solution to the 
                light time equation. For the solar system bodies, a solution is usually 
                found in 3 iterations.

        - `axescenter`: ID or instance of the center point for `axes`. This parameter is used 
                only when `axes` have an orientation that depends on time. In these cases, 
                the point is used to find the time at which the orientation of the axes 
                should be computed. It defaults to `from`.
        
        !!! note 
            If the `StellarAberrationCorrection` is applied, the frame system must be at 
            least one order higher than that of the requested transformation.

        ### See also 
        See also [`LightTime`](@ref), [`StellarAberration`](@ref) and [`vector6`](@ref).

        ### References 
        - CSPICE [Library](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkez_c.html)
        """
        @inline function ($pfun)(frames::FrameSystem{<:Any, <:Any, S}, from, to, axes, 
            ep::Epoch{S}, ltcorr::AbstractLightTimeCorrection, dir::Int; kwargs...) where {S}

            $(pfun)(frames, from, to, axes, Tempo.j2000s(ep), ltcorr, dir; kwargs...)
        end

        """
            $($pfun)(frame, from, to, axes, t::Number, ltcorr, dir; <keyword arguments>) 

        Compute a light-time corrected $(3*$order)-elements state vector of a target point 
        relative to an observing point, in a given set of axes, at the desired time `t`,  
        expressed in seconds since [`J2000`](@ref), using the aberration flag `ltcorr` and 
        the direction `dir`.
        """
        function ($pfun)(frames::FrameSystem{O, T}, from, to, axes, t::Number, 
                    ltc::AbstractLightTimeCorrection, dir::Int; 
                    iters=1, axescenter=nothing) where {O, T}

            if O < $order
                throw(ErrorException("Insufficient frame system order: "*
                    "transformation requires at least order $($order)."))
            end

            if !(dir in (-1, 1))
                throw(
                    ArgumentError(
                        "$dir is an invalid direction. Only -1 (Reception) and 1 "*
                        "(Transmission) are accepted.")
                )
            end
        
            # Retrieve aliased point ids
            idfrom = point_alias(from)
            idto = point_alias(to) 

            idfrom == idto && return @SVector zeros(T, 3*$order)
            
            # in absence of a specified axes center, the observer location is used. 
            axc = isnothing(axescenter) ? idfrom : point_alias(axescenter)

            ltp = LTProperties(axes_alias(axes), axc, dir, iters)
            $(ltcorr)(frames, ltc, ltp, idfrom, idto, t)
            
        end

    end

end


# ---------------------------------------------
# UTILITIES
# ---------------------------------------------
""" 
    update_point!(frames, point, stv::AbstractVector, epoch::Epoch)
"""
function update_point!(frames::FrameSystem{O, T, S}, point, stv::AbstractVector{T}, 
            epoch::Epoch{S}) where {O, T, S}
    update_point!(frames, point, stv, Tempo.j2000s(epoch))
end


""" 
    update_point!(frames::FrameSystem, point, stv::AbstractVector, time)

Update the state vector of `point` at the input `time` in `frames`. The only 
accepted length for the input vector `stv` are 3, 6, 9 or 12. The order is automatically 
inferred from the vector length.

### Examples 
```jldoctest
julia> FRAMES = FrameSystem{2, Float64}();
  
julia> @axes ICRF 1  

julia> add_axes_inertial!(FRAMES, ICRF)

julia> @point Origin 0

julia> @point Satellite 1 

julia> add_point_root!(FRAMES, Origin, ICRF)

julia> add_point_updatable!(FRAMES, Satellite, Origin, ICRF)

julia> y = [10000., 200., 300.]

julia> update_point!(FRAMES, Satellite, y, 0.1)

julia> vector3(FRAMES, Origin, Satellite, ICRF, 0.1)
3-element SVector{3, Float64} with indices SOneTo(3):
 10000.0
   200.0
   300.0

julia> vector3(FRAMES, Origin, Satellite, ICRF, 0.2)
ERROR: UpdatablePoint with NAIFId = 1 has not been updated at time 0.2 for order 1

julia> vector6(FRAMES, Origin, Satellite, ICRF, 0.1)
ERROR: UpdatablePoint with NAIFId = 1 has not been updated at time 0.2 for order 2
```
### See also 
See also [`add_point_updatable!`](@ref)
"""
function update_point!(frames::FrameSystem{O, T}, point, stv::AbstractVector{T}, 
            time::T) where {O, T}

    NAIFId = point_alias(point) 

    !has_point(frames, NAIFId) && throw(ErrorException(
        "point with NAIF ID $NAIFId is not contained in the frame system."))

    ne = length(stv) 
    ne > 3*O && throw(ArgumentError("state vector order greater than frame system order."))
  
    if !(ne in (3, 6, 9, 12)) 
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
@generated function _get_comp_axes_vector3(frame::FrameSystem{O, T}, v::ComputableAxesVector, 
                        axesid::Int, t::Number) where {O, T}

    expr = :(tuple($([0 for i in 1:3(O-1)]...)))

    return quote 
        @inbounds begin 
            if v.order == 1 
                stv = vector3(frame, v.from, v.to, axesid, t)        
                return SA[stv..., $(expr)...]
            elseif v.order == 2 
                stv = vector6(frame, v.from, v.to, axesid, t)
                return SA[stv[4], stv[5], stv[6], $(expr)...]
            else 
                stv = vector9(frame, v.from, v.to, axesid, t)
                return SA[stv[7], stv[8], stv[9], $(expr)...]
            end
        end
    end 
end


@generated function _get_comp_axes_vector6(frame::FrameSystem{O, T}, v::ComputableAxesVector, 
                        axesid::Int, t::Number) where {O, T}

    expr = :(tuple($([0 for i in 1:3(O-2)]...)))

    return quote 
        @inbounds begin 
            if v.order == 1 
                stv = vector6(frame, v.from, v.to, axesid, t)        
                return SA[stv..., $(expr)...]
            elseif v.order == 2 
                stv = vector9(frame, v.from, v.to, axesid, t)
                return SA[stv[4], stv[5], stv[6], stv[7], stv[8], stv[9], $(expr)...]
            else 
                stv = vector12(frame, v.from, v.to, axesid, t)
                return SA[stv[7], stv[8], stv[9], stv[10], stv[11], stv[12], $(expr)...]
            end
        end
    end 
end


@generated function _get_comp_axes_vector9(frame::FrameSystem{O, T}, v::ComputableAxesVector, 
                        axesid::Int, t::Number) where {O, T}

    expr = :(tuple($([0 for i in 1:3(O-3)]...)))

    return quote 
        @inbounds begin 
            if v.order == 1 
                stv = vector9(frame, v.from, v.to, axesid, t)        
                return SA[stv..., $(expr)...]

            elseif v.order == 2 # v.order = 2
                stv = vector12(frame, v.from, v.to, axesid, t)
                return SA[stv[4], stv[5], stv[6], stv[7], stv[8], stv[9], 
                         stv[10], stv[11], stv[12], $(expr)...]
            else 
                throw(ErrorException("unable to compute a vector of order 5 (jounce)."))
            end
        end
    end 
end


@generated function _get_comp_axes_vector12(frame::FrameSystem{O, T}, v::ComputableAxesVector, 
                        axesid::Int, t::Number) where {O, T}

    return quote 
        @inbounds begin 
            if v.order == 1 
                return vector12(frame, v.from, v.to, axesid, t)        
            elseif v.order == 2 # v.order = 2
                throw(ErrorException("unable to compute a vector of order 5 (snap)."))

            else 
                throw(ErrorException("unable to compute a vector of order 6 (crackle)."))
            end
        end
    end 

end

