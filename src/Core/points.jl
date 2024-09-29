""" 
    add_point!(frames, name, id, axesid, class, funs, parentid=nothing)

Create and add a new point node `name` to `frames` based on the input parameters. 

### Inputs 
- `frames` -- Target frame system 
- `name` -- Point name, must be unique within `frames` 
- `id` -- Point ID, must be unique within `frames`
- `axes` -- ID/Name of the axes in which the state vector of the point is expressed. 
- `funs` -- `FramePointFunctions` object storing the functions to update the state 
            vectors of the point.
- `parentid` -- NAIF ID of the parent point. Not required only for the root point.

!!! warning 
    This is a low-level function and is NOT meant to be directly used. Instead, to add a point 
    to the frame system, see [`add_point_dynamical!`](@ref), [`add_point_fixedoffset!`](@ref)
    and [`add_point_root!`](@ref).
"""
function add_point!(
    frames::FrameSystem{O,T}, name::Symbol, id::Int, axes,
    funs::FramePointFunctions{O,T}=FramePointFunctions{O,T}(), parentid=nothing
) where {O,T<:Number}

    if has_point(frames, id)
        # Check point with the same id already registered 
        throw(
            ArgumentError(
                "A point with ID $id is already registered in the input frame system.",
            ),
        )
    end

    # Check point with the same name does not already exist 
    if name in map(x -> x.name, points_graph(frames).nodes)
        throw(
            ArgumentError(
                "A point with name=$name is already registed in the input frame system"
            ),
        )
    end

    # Check if the given axes are known in the FrameSystem
    axesid = axes_id(frames, axes)
    if !has_axes(frames, axesid)
        throw(
            ArgumentError(
                "Axes with ID $axesid are not registered in the input frame system"
            ),
        )
    end

    if isnothing(parentid)
        # If a root-point exists, check that a parent has been specified 
        if !isempty(points_graph(frames))
            throw(
                ArgumentError(
                    "A parent point is required because the input frame system " *
                    "already contains a root-point.",
                ),
            )
        end

        parentid = id # Root-point has parentid = id

    else
        # Check that the parent point is registered in frames 
        if !has_point(frames, parentid)
            throw(
                ArgumentError(
                    "The specified parent point with id $parentid is not " *
                    "registered in the input frame system.",
                ),
            )
        end
    end

    # Creates point node 
    pnt = FramePointNode{O,T}(name, id, parentid, axesid, funs)

    # Insert new point in the graph
    add_point!(frames, pnt)

    # Connect the new point to the parent point in the graph 
    !isnothing(parentid) && add_edge!(points_graph(frames), parentid, id)

    return nothing
end

"""
    add_point_fixedoffset!(frames, name, id, parent, axes, offset::AbstractVector)

Add `point` as a fixed-offset point to `frames`. 

Fixed points are those whose positions have a constant `offset` with respect their `parent` 
points in the given set of `axes`. Thus, points eligible for this class must have null 
velocity and acceleration with respect to `parent`.
"""
function add_point_fixedoffset!(
    frames::FrameSystem{O,T}, name::Symbol, id::Int, parent, ax,
    offset::AbstractVector{N}
) where {O,N,T}

    if length(offset) != 3
        throw(
            DimensionMismatch(
                "The offset vector should have length 3, but has $(length(offset))."
            ),
        )
    end

    tr = Translation{O}(SVector(offset...))
    funs = FramePointFunctions{O,T}(t -> tr)

    return add_point!(
        frames, name, id, axes_id(frames, ax), funs, point_id(frames, parent)
    )
end

""" 
    add_point_dynamical!(frames, name, id, parent, axes, fun, δfun=nothing, δ²fun=nothing, δ³fun=nothing)

Add `point` as a time point to `frames`. The state vector for these points depends only on 
time and is computed through the custom functions provided by the user. 

The input functions must accept only time as argument and their outputs must be as follows: 

- **fun**: return a 3-elements vector: position
- **δfun**: return a 6-elements vector: position and velocity
- **δ²fun**: return a 9-elements vector: position, velocity and acceleration
- **δ³fun**: return a 12-elements vector: position, velocity, acceleration and jerk

If `δfun`, `δ²fun` or `δ³fun` are not provided, they are computed with automatic differentiation. 

!!! warning 
    It is expected that the input functions and their ouputs have the correct signature. This 
    function does not perform any checks on whether the returned vectors have the appropriate 
    dimensions. 
"""
function add_point_dynamical!(
    frames::FrameSystem{O,T}, name::Symbol, id::Int, parent, ax, fun,
    δfun=nothing, δ²fun=nothing, δ³fun=nothing
) where {O,T}

    for (order, fcn) in enumerate([δfun, δ²fun, δ³fun])
        if (O < order + 1 && !isnothing(fcn))
            @warn "ignoring $fcn, frame system order is less than $(order+1)"
        end
    end

    funs = FramePointFunctions{O,T}(
        t -> Translation{O}(fun(t)),

        # First derivative
        if isnothing(δfun)
            t -> Translation{O}(vcat(fun(t), D¹(fun, t)))
        else
            t -> Translation{O}(δfun(t))
        end,

        # Second derivative
        if isnothing(δ²fun)
            (
                if isnothing(δfun)
                    t -> Translation{O}(vcat(fun(t), D¹(fun, t), D²(fun, t)))
                else
                    t -> Translation{O}(vcat(δfun(t), D²(fun, t)))
                end
            )
        else
            t -> Translation{O}(δ²fun(t))
        end,

        # Third derivative 
        if isnothing(δ³fun)
            (
                if isnothing(δ²fun)
                    (
                        if isnothing(δfun)
                            t -> Translation{O}(vcat(fun(t), D¹(fun, t), D²(fun, t), D³(fun, t)))
                        else
                            t -> Translation{O}(vcat(δfun(t), D²(fun, t), D³(fun, t)))
                        end
                    )
                else
                    t -> Translation{O}(vcat(δ²fun(t), D³(fun, t)))
                end
            )
        else
            t -> Translation{O}(δ³fun(t))
        end,
    )

    return add_point!(
        frames, name, id, axes_id(frames, ax), funs, point_id(frames, parent)
    )
end