for (order, axfun, _axfun, pfun, _pfun, _pfwd, _pbwd) in zip(
    (1, 2, 3, 4),
    (:rotation3, :rotation6, :rotation9, :rotation12),
    (:_rotation3, :_rotation6, :_rotation9, :_rotation12),
    (:vector3, :vector6, :vector9, :vector12),
    (:_vector3, :_vector6, :_vector9, :_vector12),
    (:_vector3_forward, :_vector6_forward, :_vector9_forward, :_vector12_forward),
    (:_vector3_backward, :_vector6_backward, :_vector9_backward, :_vector12_backward),
)

    # --------------------------------------------------------------------------------------
    # Axes transformations
    # --------------------------------------------------------------------------------------

    @eval begin 

        @inline function ($axfun)(
            ::FrameSystem{<:Any, <:Any, S1}, from, to, ::Epoch{S2}
        ) where {S1, S2}
            throw(ArgumentError("Incompatible epoch timescale: expected $S1, found $S2."))
        end

        @inline function ($axfun)(
            frame::FrameSystem{<:Any,<:Any, S}, from, to, ep::Epoch{S}
        ) where {S}
            return $(axfun)(frame, from, to, Tempo.j2000s(ep))
        end

        function ($axfun)(frame::FrameSystem{O, T}, from::Int, to::Int, t::Number) where {O,T}
            if O < $order
                throw(
                    ErrorException(
                        "insufficient frame system order: " *
                        "transformation requires at least order $($order).",
                    ),
                )
            end

            from == to && return Rotation{$order}(T(1)I)

            # # Check to ensure that the two axes are stored in the frame system
            # for id in (from, to)
            #     if !has_axes(frame, id)
            #         throw(
            #             ErrorException(
            #                 "axes with ID $id are not registered in the frame system."
            #             )
            #         )
            #     end 
            # end
            return $(_axfun)(frame, get_path(get_axes(frame), from, to), t)
        end

        # Low-level function to parse a path of axes and chain their rotations 
        @inbounds function ($_axfun)(frame::FrameSystem, path::Vector{Int}, t::Number)
            f1 = get_mappednode(get_axes(frame), path[1])
            f2 = get_mappednode(get_axes(frame), path[2])
            rot = $(_axfun)(f1, f2, t)

            for i in 2:(length(path) - 1)
                f1 = f2
                f2 = get_mappednode(get_axes(frame), path[i + 1])
                rot = $(_axfun)(f1, f2, t) * rot
            end
            return rot
        end

        # Low-level function to compute the rotation between two axes
        @inline function ($_axfun)(from::FrameAxesNode, to::FrameAxesNode, t::Number)
            return if from.id == to.parentid
                return $(_axfun)(to, t)
            else
                return inv($(_axfun)(from, t))
            end
        end

        @inline function($_axfun)(ax::FrameAxesNode, t::Number)
            R = ax.f[$order](t) 
            return Rotation{$order}(R)
        end

    end

    # --------------------------------------------------------------------------------------
    # Points transformations
    # --------------------------------------------------------------------------------------

    @eval begin 

        @inline function ($pfun)(
            ::FrameSystem{<:Any,<:Any,S1}, from, to, axes, ::Epoch{S2}
        ) where {S1,S2}
            throw(ArgumentError("Incompatible epoch timescale: expected $S1, found $S2."))
        end

        @inline function ($pfun)(
            frame::FrameSystem{<:Any,<:Any,S}, from, to, axes, ep::Epoch{S}
        ) where {S}
            return $(pfun)(frame, from, to, axes, Tempo.j2000s(ep))
        end

        function ($pfun)(frame::FrameSystem{O, N}, from::Int, to::Int, axes::Int, t::Number) where {O, N}
            if O < $order
                throw(
                    ErrorException(
                        "insufficient frame system order: " *
                        "transformation requires at least order $($order).",
                    ),
                )
            end

            from == to && return @SVector zeros(N, 3 * $order)

            # Check to ensure that the two points are registerd
            for id in (from, to)
                if !has_point(frame, id)
                    throw(
                        ErrorException("point with ID $id is not registered in the frame system.")
                    )
                end 
            end

            # Check that the ouput axes are registered 
            if !has_axes(frame, axes)
                throw(
                    ErrorException("axes with ID $axes are not registered in the frame system.")
                )
            end

            return $(_pfun)(frame, get_path(get_points(frame), from, to), axes, t)
        end

        function ($_pfun)(frame::FrameSystem, path::Vector{Int}, axes::Int, t::Number)
            @inbounds p1 = get_mappednode(get_points(frame), path[1])
            @inbounds p2 = get_mappednode(get_points(frame), path[end])

            if length(path) == 2 
                # This handles all the cases where you don't need to chain any transformations
                axid, stv = ($_pfun)(p1, p2, t)
                if axid != axes
                    return $(axfun)(frame, axid, axes, t) * stv 
                end
                return stv
            elseif axes == p1.axesid 
                # backward pass 
                return $(_pbwd)(frame, p2, path, t)
            elseif axes == p2.axesid 
                # forward pass 
                return $(_pfwd)(frame, p1, path, t)
            else 
                # Optimising this transformation would probably demand a significant 
                # portion of time with respect to the time required by the whole transformation
                # therefore forward pass is used without any optimisation
                return $(axfun)(frame, p2.axesid, axes, t) * $(_pfwd)(frame, path, p1, t)
            end
        end

        @inbounds function ($_pfwd)(frame::FrameSystem, p1::FramePointNode, path::Vector{Int}, t::Number)
            p2 = get_mappednode(get_points(frame), path[2])
            axid, stv = ($_pfun)(p1, p2, t)
            for i in 2:(length(path)-1)
                p1 = p2 
                p2 = get_mappednode(get_points(frame), path[i+1])
                ax2id, stv2 = ($_pfun)(p1, p2, t)

                # Rotates previous vector to p2's axes
                if ax2id != axid 
                    stv = ($axfun)(frame, axid, ax2id, t) * stv
                end

                axid = ax2id 
                stv += stv2
            end
            return stv 
        end

        @inbounds function ($_pbwd)(frame::FrameSystem, p1::FramePointNode, path::Vector{Int}, t::Number)
            p2 = get_mappednode(get_points(frame), path[end-1])
            axid, stv = ($_pfun)(p1, p2, t)
            for i in 2:(length(path)-1)
                p1 = p2 
                p2 = get_mappednode(get_points(frame), path[end-i])
                ax2id, stv2 = ($_pfun)(p1, p2, t)

                # Rotates previous vector to p2's axes
                if ax2id != axid 
                    stv = ($axfun)(frame, axid, ax2id, t) * stv
                end

                axid = ax2id 
                stv += stv2
            end
            return -stv 
        end

        @inbounds function ($_pfun)(from::FramePointNode, to::FramePointNode, t::Number)
            if from.id == to.parentid 
                return to.axesid, $(_pfun)(to, t)
            else 
                return from.axesid, -$(_pfun)(from, t)
            end
        end

        @inbounds function ($_pfun)(p::FramePointNode, t::Number)
            stv = p.f[$order](t)
            D = 3 * $order
            @views return SA[(stv[1:D])...]
        end

    end

end