
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
    frames::FrameSystem{O,N}, name::Symbol, axes, fun::Function,
    δfun=nothing, δ²fun=nothing, δ³fun=nothing
) where {O,N}

    for (order, fcn) in enumerate([δfun, δ²fun, δ³fun])
        if (O < order + 1 && !isnothing(fcn))
            @warn "ignoring $fcn, frame system order is less than $(order+1)"
        end
    end

    funs = DirectionFunctions{O,N}(
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
                            t -> Direction{O}(vcat(fun(t), D¹(fun, t), D²(fun, t), D³(fun, t)))
                        else
                            t -> Direction{O}(vcat(δfun(t), D²(fun, t), D³(fun, t)))
                        end
                    )
                else
                    t -> Direction{O}(vcat(δ²fun(t), D³(fun, t)))
                end
            )
        else
            t -> Translation{O}(δ³fun(t))
        end,
    )

    axid = axes_id(frames, axes)
    dir = DirectionDefinition{O,N}(name, length(directions(frames)) + 1, axid, funs)
    push!(directions(frames), Pair(name, dir))
    nothing
end
