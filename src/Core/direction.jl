
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

function add_direction!(
    frames::FrameSystem{O, N}, name::Symbol,  
    fun::Function, δfun = nothing, δ²fun = nothing, δ³fun = nothing,
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
    return add_direction!(frames, name, funs)
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

function add_direction_velocity!(
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

    fun  = t->SVectorNT{3O, N}(@views(vector6(frames, from, to, ax, t)[4:end]))
    dfun = t->SVectorNT{3O, N}(@views(vector9(frames, from, to, ax, t)[4:end]))
    ddfun= t->SVectorNT{3O, N}(@views(vector12(frames, from, to, ax, t)[4:end]))

    return add_direction!(frames, name, fun, dfun , ddfun)
end

function add_direction_orthogonal!(
    frames::FrameSystem{O, N}, name::Symbol, dir1::Symbol, dir2::Symbol
) where {O, N}

    for d in (dir1, dir2)
        if !(d in directions(frames))
            throw(
                ErrorException(
                    "no direction with name $d registered in the given frame system"
                )
            )
        end
    end

    fun    = t->SVectorNT{3O, N}(
        cross3(direction3(frames, dir1, t), direction3(frames, dir2, t))
    )
    dfun   = t->SVectorNT{3O, N}(
        cross6(direction6(frames, dir1, t), direction6(frames, dir2, t))
    ) 
    ddfun  = t->SVectorNT{3O, N}(
        cross9(direction9(frames, dir1, t), direction9(frames, dir2, t))
    ) 
    dddfun = t->SVectorNT{3O, N}(
        cross12(direction12(frames, dir1, t), direction12(frames, dir2, t))
    ) 

    return add_direction!(frames, name, fun, dfun, ddfun, dddfun)
end

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
    fun    = t->SVectorNT{3O, N}(unitvec(direction3(frames, dir, t)))
    dfun   = t->SVectorNT{3O, N}(_normalize6(frames, dir, t)) 
    ddfun  = t->SVectorNT{3O, N}(_normalize9(frames, dir, t))
    dddfun = t->SVectorNT{3O, N}(_normalize12(frames, dir, t))
    return add_direction!(frames, name, fun, dfun, ddfun, dddfun)
end

function _normalize6(frames, dir, t)
    d = direction6(frames, dir, t)
    @views nd = vcat(unitvec(d[1:3]), δunitvec(d[4:6]))
    return nd
end

function _normalize9(frames, dir, t)
    d = direction9(frames, dir, t)
    @views nd = vcat(unitvec(d[1:3]), δunitvec(d[4:6]), δ²unitvec(d[7:9]))
    return nd
end

function _normalize12(frames, dir, t)
    d = direction12(frames, dir, t)
    @views nd = vcat(unitvec(d[1:3]), δunitvec(d[4:6]), δ²unitvec(d[7:9]), δ³unitvec(d[10:12]))
    return nd
end