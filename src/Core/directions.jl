
"""
    add_direction!(frames, name::Symbol, axesid, funs)

Add a new direction node to `frames`.

### Inputs 
- `frames` -- Target frame system 
- `name` -- Direction name, must be unique within `frames` 
- `axesid` -- ID of the axes the direction is expressed in
- `funs` -- `DirectionFunctions` object storing the functions to compute the direction and, 
            eventually, its time derivatives. It must match the type and order of `frames`.
"""
function add_direction!(
    frames::FrameSystem{O, N}, name::Symbol, axesid::Int, funs::DirectionFunctions{O, N}
) where {O, N <: Number}
    if name in directions(frames)
        throw(
            ArgumentError(
                "A direction with name $name is already registered in the input frame system.",
            ),
        )
    end

    dir = Direction{O, N}(name, length(directions(frames))+1, axesid, funs)
    push!(directions_map(frames), Pair(name, dir))
    nothing
end

"""
    add_direction!(frames, name::Symbol, axes, fun, δfun=nothing, δ²fun=nothing, δ³fun=nothing)

Add a new direction node to `frames`. The orientation of these direction depends only 
on time and is computed through the custom functions provided by the user. 

The input functions must accept only time as argument and their outputs must be as follows: 

- `fun`: return a direction vector.
- `δfun`: return a direction vector and its 1st order time derivative.
- `δ²fun`: return a direction vector and its 1st and 2nd order time derivatives.
- `δ³fun`: return a direction vector and its 1st, 2nd and 3rd order time derivatives.

If `δfun`, `δ²fun` or `δ³fun` are not provided, they are computed via automatic differentiation.

!!! warning 
    It is expected that the input functions and their outputs have the correct signature. This 
    function does not perform any checks on the output types. 
"""
function add_direction!(
    frames::FrameSystem{O, N}, name::Symbol, ax, fun::Function, 
    δfun = nothing, δ²fun = nothing, δ³fun = nothing
) where {O, N}

    for (order, fcn) in enumerate([δfun, δ²fun, δ³fun])
        if (O < order + 1 && !isnothing(fcn))
            @warn "ignoring $fcn, frame system order is less than $(order+1)"
        end
    end

    funs = DirectionFunctions{O, N}(
        t -> SVectorNT{3O, N}(fun(t)),

        # First derivative
        if isnothing(δfun)
            t -> SVectorNT{3O, N}(vcat(fun(t), D¹(fun, t)))
        else
            t -> SVectorNT{3O, N}(δfun(t))
        end,

        # Second derivative
        if isnothing(δ²fun)
            (
                if isnothing(δfun)
                    t -> SVectorNT{3O, N}(vcat(fun(t), D¹(fun, t), D²(fun, t)))
                else
                    t -> SVectorNT{3O, N}(vcat(δfun(t), D²(fun, t)))
                end
            )
        else
            t -> SVectorNT{3O, N}(δ²fun(t))
        end,

        # Third derivative 
        if isnothing(δ³fun)
            (
                if isnothing(δ²fun)
                    (
                        if isnothing(δfun)
                            t -> SVectorNT{3O, N}(vcat(fun(t), D¹(fun, t), D²(fun, t), D³(fun, t)))
                        else
                            t -> SVectorNT{3O, N}(vcat(δfun(t), D²(fun, t), D³(fun, t)))
                        end
                    )
                else
                    t -> SVectorNT{3O, N}(vcat(δ²fun(t), D³(fun, t)))
                end
            )
        else
            t -> SVectorNT{3O, N}(δ³fun(t))
        end,
    )
    return add_direction!(frames, name, axes_id(frames, ax), funs)
end

"""
    add_direction_fixed!(frames, name, axes, offset::AbstractVector)

Add a fixed direction to `frames`.
"""
function add_direction_fixed!(
    frames::FrameSystem{O, N}, name::Symbol, ax, offset::AbstractVector{T}
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

    return add_direction!(frames, name, ax, funs)
end

"""
    add_direction_position!(frames, name::Symbol, origin, target, axes)

Add a direction based on the position vector from `origin` to `target` in the specified `axes`.
"""
function add_direction_position!(
    frames::FrameSystem{O, N}, name::Symbol, from, to, ax
) where {O, N}

    fromid, toid = point_id(frames, from),  point_id(frames, to)
    axid = axes_id(frames, ax)

    for id in (fromid, toid)
        !has_point(frames, id) && throw(
            ArgumentError(
                "no points with id $id registered in the given frame system."
            )
        )
    end

    !has_axes(frames, axid) && throw(
        ArgumentError(
            "no axes with id $ax registered in the given frame system"
        )
    )

    funs = DirectionFunctions{O, N}(
        t->SVectorNT{3O, N}(vector3(frames, fromid, toid, axid, t)),
        t->SVectorNT{3O, N}(vector6(frames, fromid, toid, axid, t)),
        t->SVectorNT{3O, N}(vector9(frames, fromid, toid, axid, t)),
        t->SVectorNT{3O, N}(vector12(frames, fromid, toid, axid, t))
    )
    return add_direction!(frames, name, axid, funs)
end

"""
    add_direction_velocity!(frames, name::Symbol, origin, target, axes)

Add a direction based on the velocity vector from `origin` to `target` in the specified `axes`.
"""
function add_direction_velocity!(
    frames::FrameSystem{O, N}, name::Symbol, from, to, ax
) where {O, N}

    fromid, toid = point_id(frames, from),  point_id(frames, to)
    axid = axes_id(frames, ax)

    for id in (fromid, toid)
        !has_point(frames, id) && throw(
            ArgumentError(
                "no points with id $id registered in the given frame system."
            )
        )
    end

    !has_axes(frames, axid) && throw(
        ArgumentError(
            "no axes with id $ax registered in the given frame system"
        )
    )

    fun  = t->SVectorNT{3O, N}(@views(vector6(frames, fromid, toid, axid, t)[4:end]))
    dfun = t->SVectorNT{3O, N}(@views(vector9(frames, fromid, toid, axid, t)[4:end]))
    ddfun= t->SVectorNT{3O, N}(@views(vector12(frames, fromid, toid, axid, t)[4:end]))

    return add_direction!(frames, name, axid, fun, dfun , ddfun)
end

"""
    add_direction_orthogonal!(frames, name::Symbol, dir1, dir2)

Add a direction as the cross product between two existing directions (i.e. `dir1` and `dir2`).
"""
function add_direction_orthogonal!(
    frames::FrameSystem{O, N}, name::Symbol, dir1::Symbol, dir2::Symbol, ax
) where {O, N}

    for d in (dir1, dir2)
        if !(d in directions(frames))
            throw(
                ArgumentError(
                    "no direction with name $d registered in the given frame system"
                )
            )
        end
    end

    axid = axes_id(frames, ax)

    !has_axes(frames, axid) && throw(
        ArgumentError(
            "no axes with id $ax registered in the given frame system"
        )
    )

    fun    = t->SVectorNT{3O, N}(
        cross3(direction3(frames, dir1, axid, t), direction3(frames, dir2, axid, t))
    )
    dfun   = t->SVectorNT{3O, N}(
        cross6(direction6(frames, dir1, axid, t), direction6(frames, dir2, axid, t))
    ) 
    ddfun  = t->SVectorNT{3O, N}(
        cross9(direction9(frames, dir1, axid, t), direction9(frames, dir2, axid, t))
    ) 
    dddfun = t->SVectorNT{3O, N}(
        cross12(direction12(frames, dir1, axid, t), direction12(frames, dir2, axid, t))
    ) 

    return add_direction!(frames, name, axid, fun, dfun, ddfun, dddfun)
end

"""
    add_direction_normalize!(frames, name::Symbol, dir)

Add a direction as the normalized version of `dir`.
"""
function add_direction_normalize!(
    frames::FrameSystem{O, N}, name::Symbol, dir::Symbol
) where {O, N}
    if !(dir in directions(frames))
        throw(
            ErrorException(
                "no direction with name $dir registered in the given frame system"
            )
        )
    end
    axid = directions_map(frames)[dir].axesid

    fun    = t->SVectorNT{3O, N}(unitvec(direction3(frames, dir, axid, t)))
    dfun   = t->SVectorNT{3O, N}(_normalize6(frames, dir, axid, t)) 
    ddfun  = t->SVectorNT{3O, N}(_normalize9(frames, dir, axid, t))
    dddfun = t->SVectorNT{3O, N}(_normalize12(frames, dir, axid, t))
    return add_direction!(frames, name, fun, dfun, ddfun, dddfun)
end

function _normalize6(frames, dir, axid, t)
    d = direction6(frames, dir, axid, t)
    @views nd = vcat(unitvec(d[1:3]), δunitvec(d[4:6]))
    return nd
end

function _normalize9(frames, dir, axid, t)
    d = direction9(frames, dir, axid, t)
    @views nd = vcat(unitvec(d[1:3]), δunitvec(d[4:6]), δ²unitvec(d[7:9]))
    return nd
end

function _normalize12(frames, dir, axid, t)
    d = direction12(frames, dir, axid, t)
    @views nd = vcat(unitvec(d[1:3]), δunitvec(d[4:6]), δ²unitvec(d[7:9]), δ³unitvec(d[10:12]))
    return nd
end