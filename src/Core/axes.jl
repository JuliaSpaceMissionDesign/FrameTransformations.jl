"""
    add_axes!(frames, name::Symbol, id::Int, class::Int, funs, parentid)

Add a new axes node to `frames`.

### Inputs 
- `frames` -- Target frame system 
- `name` -- Axes name, must be unique within `frames` 
- `id` -- Axes ID, must be unique within `frames`
- `funs` -- `FrameAxesFunctions` object storing the functions to compute the DCM and, 
            eventually, its time derivatives. 
- `parentid` -- Axes ID of the parent axes. Not required only for the root axes.

!!! warning 
    This is a low-level function and is NOT meant to be directly used. Instead, to add a set of
    axes to the frame system, see [`add_axes_projected!`](@ref), [`add_axes_rotating!`](@ref) 
    and [`add_axes_fixedoffset!`](@ref).
"""
function add_axes!(
    frames::FrameSystem{O,T}, name::Symbol, id::Int,
    funs::FrameAxesFunctions{O,T}=FrameAxesFunctions{O,T}(),
    parentid=nothing
) where {O,T<:Number}

    if has_axes(frames, id)
        # Check if a set of axes with the same ID is already registered within 
        # the given frame system 
        throw(
            ArgumentError(
                "Axes with ID $id are already registered in the frame system."
            ),
        )
    end

    if name in map(x -> x.name, axes_graph(frames).nodes)
        # Check if axes with the same name also do not already exist
        throw(
            ArgumentError(
                "Axes with name=$name are already registered in the frame system."
            ),
        )
    end

    if !isnothing(parentid)
        # Check if the root axes is not present
        isempty(axes_graph(frames)) && throw(ArgumentError("Missing root axes."))

        # Check if the parent axes are registered in frame 
        if !has_axes(frames, parentid)
            throw(
                ArgumentError(
                    "The specified parent axes with ID $parentid are not " *
                    "registered in the frame system.",
                ),
            )
        end
    else
        # Check if axes are already present 
        !isempty(axes_graph(frames)) && throw(ArgumentError("Root axes already registed."))

        # Root axes
        parentid = id
    end


    # Create point
    node = FrameAxesNode{O,T}(name, id, parentid, funs)

    # Insert new point in the graph
    add_axes!(frames, node)

    # Connect the new axes to the parent axes in the graph 
    !isnothing(parentid) && add_edge!(axes_graph(frames), parentid, id)

    return nothing
end

"""
    add_axes_fixedoffset!(frames, name::Symbol, id::Int, parent, dcm:DCM)
   
Add axes `name` with id `id` to `frames` with a fixed-offset from `parent`. 
Fixed offset axes have a constant orientation with respect to their `parent` axes, 
represented by `dcm`, a Direction Cosine Matrix (DCM).

### See also 
See also [`add_axes!`](@ref).
"""
function add_axes_fixedoffset!(
    frames::FrameSystem{O,T}, name::Symbol, id::Int, parent, dcm::DCM{T}
) where {O,T}

    funs = FrameAxesFunctions{O,T}(t -> Rotation{O}(dcm))
    add_axes!(frames, name, id, funs, axes_id(frames, parent))
end

"""
    add_axes_projected!(frames, name, id, parent, fun)

Add inertial axes `name` and id `id` as a set of projected axes to `frames`. The axes relation 
to the `parent` axes are given by a `fun`. 

Projected axes are similar to rotating axis, except that all the positions, velocity, etc ... 
are rotated by the 0-order rotation (i.e. the derivatives of the rotation matrix are null, 
despite the rotation depends on time).

### See also 
See also [`add_axes!`](@ref).
"""
function add_axes_projected!(
    frames::FrameSystem{O,T}, name::Symbol, id::Int, parent, fun::Function
) where {O,T}
    funs = FrameAxesFunctions{O,T}(t -> Rotation{O}(fun(t)))
    add_axes!(frames, name, id, funs, axes_id(frames, parent))
end

"""
    add_axes_rotating!(frames, name::Symbol, id::Int, parent, fun, δfun=nothing, 
        δ²fun=nothing, δ³fun=nothing)
   
Add `axes` as a set of rotating axes to `frames`. The orientation of these axes depends only 
on time and is computed through the custom functions provided by the user. 

The input functions must accept only time as argument and their outputs must be as follows: 

- `fun`: return a Direction Cosine Matrix (DCM).
- `δfun`: return the DCM and its 1st order time derivative.
- `δ²fun`: return the DCM and its 1st and 2nd order time derivatives.
- `δ³fun`: return the DCM and its 1st, 2nd and 3rd order time derivatives.

If `δfun`, `δ²fun` or `δ³fun` are not provided, they are computed via automatic differentiation.

!!! warning 
    It is expected that the input functions and their outputs have the correct signature. This 
    function does not perform any checks on the output types. 
"""
function add_axes_rotating!(
    frames::FrameSystem{O,T}, name::Symbol, id::Int, parent, fun,
    δfun=nothing, δ²fun=nothing, δ³fun=nothing,
) where {O,T}

    for (order, fcn) in enumerate([δfun, δ²fun, δ³fun])
        if (O < order + 1 && !isnothing(fcn))
            @warn "ignoring $fcn, frame system order is less than $(order+1)"
        end
    end

    funs = FrameAxesFunctions{O,T}(
        t -> Rotation{O}(fun(t)),

        # First derivative 
        if isnothing(δfun)
            t -> Rotation{O}(fun(t), D¹(fun, t))
        else
            t -> Rotation{O}(δfun(t))
        end,

        # Second derivative 
        if isnothing(δ²fun)
            (
                if isnothing(δfun)
                    t -> Rotation{O}(fun(t), D¹(fun, t), D²(fun, t))
                else
                    t -> Rotation{O}(δfun(t)..., D²(fun, t))
                end
            )
        else
            t -> Rotation{O}(δ²fun(t))
        end,

        # Third derivative 
        if isnothing(δ³fun)
            (
                if isnothing(δ²fun)
                    (
                        if isnothing(δfun)
                            t -> Rotation{O}(fun(t), D¹(fun, t), D²(fun, t), D³(fun, t))
                        else
                            t -> Rotation{O}(δfun(t)..., D²(δfun, t)...)
                        end
                    )
                else
                    t -> Rotation{O}(δ²fun(t)..., D³(fun, t))
                end
            )
        else
            t -> Rotation{O}(δ³fun(t))
        end,
    )

    return add_axes!(frames, name, id, funs, axes_id(frames, parent))
end

"""
    add_axes_alias!(frames, target, alias::Symbol)

Add a name `alias` to a `target` axes registered in `frames`.
"""
function add_axes_alias!(frames::FrameSystem{O,T}, target, alias::Symbol) where {O,T}
    if !has_axes(frames, target)
        throw(
            ErrorException(
                "no axes with ID $target registered in the given frame system"
            )
        )
    end

    if alias in keys(axes_alias(frames))
        throw(
            ErrorException(
                "axes with name $alias already present in the given frame system"
            )
        )
    end

    push!(axes_alias(frames), Pair(alias, axes_id(frames, target)))
    nothing
end