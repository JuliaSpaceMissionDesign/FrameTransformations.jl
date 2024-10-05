
"""
    add_direction_position!(frames, name::Symbol, origin, target, axes)

Add a direction based on the position vector from `origin` to `target` in the specified `axes`.
"""
function add_direction_position!(
    frames::FrameSystem{O,T}, name::Symbol, from, to, axes
) where {O,T}

    fromid, toid = point_id(frames, from), point_id(frames, to)
    axid = axes_id(frames, axes)

    for id in (fromid, toid)
        !has_point(frames, id) && throw(
            ArgumentError(
                "no points with id $id registered in the given frame system."
            )
        )
    end

    !has_axes(frames, axid) && throw(
        ArgumentError(
            "no axes with id $axes registered in the given frame system"
        )
    )

    fun = t -> unitvec(vector3(frames, fromid, toid, axid, t))
    return add_direction!(frames, name, axid, fun)
end

"""
    add_direction_velocity!(frames, name::Symbol, origin, target, axes)

Add a direction based on the velocity vector from `origin` to `target` in the specified `axes`.
"""
function add_direction_velocity!(
    frames::FrameSystem{O,T}, name::Symbol, from, to, axes
) where {O,T}

    fid, tid = point_id(frames, from), point_id(frames, to)
    axid = axes_id(frames, axes)

    for id in (fid, tid)
        !has_point(frames, id) && throw(
            ArgumentError(
                "no points with id $id registered in the given frame system."
            )
        )
    end

    !has_axes(frames, axid) && throw(
        ArgumentError(
            "no axes with name $axes registered in the given frame system"
        )
    )

    fun = t -> @views(unitvec(vector6(frames, fid, tid, axid, t)[4:end]))
    return add_direction!(frames, name, axid, fun)
end

"""
    add_direction_orthogonal!(frames, name::Symbol, dir1, dir2)

Add a direction as the cross product between two existing directions (i.e. `dir1` and `dir2`).
"""
function add_direction_orthogonal!(
    frames::FrameSystem{O,T}, name::Symbol, dir1::Symbol, dir2::Symbol, axes
) where {O,T}

    for d in (dir1, dir2)
        if !(d in keys(directions(frames)))
            throw(
                ArgumentError(
                    "no direction with name $d registered in the given frame system"
                )
            )
        end
    end

    axid = axes_id(frames, axes)
    !has_axes(frames, axid) && throw(
        ArgumentError(
            "no axes with name $axes registered in the given frame system"
        )
    )

    fun = t -> unitvec(cross3(direction3(frames, dir1, axid, t), direction3(frames, dir2, axid, t)))
    return add_direction!(frames, name, axid, fun)
end

"""
    add_direction_fixed!(frames, name, axes, offset::AbstractVector)

Add a fixed direction to `frames`.
"""
function add_direction_fixed!(
    frames::FrameSystem{O,N}, name::Symbol, axes, offset::AbstractVector{T}
) where {O,N,T}

    if length(offset) != 3
        throw(
            DimensionMismatch(
                "The offset vector should have length 3, but has $(length(offset))."
            ),
        )
    end

    voffset = SVector{3,N}(offset)
    fun = t -> voffset
    return add_direction!(frames, name, axes, fun)
end