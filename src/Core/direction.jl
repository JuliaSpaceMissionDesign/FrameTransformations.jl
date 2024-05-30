
function add_direction!(
    frames::FrameSystem{O, N}, name::Symbol, funs::DirectionFunctions{O, N}
) where {O, N <: Number}

    if name in directions(frames)
        throw(
            ArgumentError(
                "A direction with name $name is already registered in the input frame system.",
            ),
        )
    end

    dir = Direction{O, N}(name, length(directions(frames))+1, funs)
    push!(get_directions(frames), Pair(name, dir))

    nothing

end

function add_direction_fixed!(
    frames::FrameSystem{O, N}, name::Symbol, offset::AbstractVector{T}
) where {O, N, T}

    if length(offset) != 3
        throw(
            DimensionMismatch(
                "The offset vector should have length 3, but has $(length(offset))."
            ),
        )
    end

    voffset = SVectorNT{3O, N}(SVector(offset...))
    funs = DirectionFunctions{O, N}(t -> voffset, t -> voffset, t -> voffset, t -> voffset)

    return add_direction!(frames, name, funs)
end

function add_direction_position!(
    frames::FrameSystem{O, N}, name::Symbol, from::Symbol, to::Symbol, ax::Symbol
) where {O, N}
    return add_direction_position!(
        frames, name, points(frames)[from], points(frames)[to], axes(frames)[ax]
    )
end

function add_direction_position!(
    frames::FrameSystem{O, N}, name::Symbol, from::Int, to::Int, ax::Int
) where {O, N}

    for id in (from, to)
        !has_point(frames, id) && throw(
            ErrorException(
                "no points with id $id registered in the given frame system."
            )
        )
    end

    !has_axes(frames, ax) && throw(
        ErrorException(
            "no axes with id $ax registered in the given frame system"
        )
    )



    funs = DirectionFunctions{O, N}(
        t->SVectorNT{3O, N}(vector3(frames, from, to, ax, t)),
        t->SVectorNT{3O, N}(vector6(frames, from, to, ax, t)),
        t->SVectorNT{3O, N}(vector9(frames, from, to, ax, t)),
        t->SVectorNT{3O, N}(vector12(frames, from, to, ax, t))
    )
    return add_direction!(frames, name, funs)
end