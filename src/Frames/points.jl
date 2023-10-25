export @point,
    add_point_root!,
    add_point_ephemeris!,
    add_point_fixed!,
    add_point_updatable!,
    add_point_dynamical!,
    point_alias

"""
    point_name(point::AbstractFramePoint)

Return the name of `point`.
"""
function point_name end

""" 
    point_id(point::AbstractFramePoint)

Return the NAIF ID associated to `point`.
"""
function point_id end

"""
    point_alias(ax::AbstractFramePoint)

Return the NAIF ID associated to the input point. 

"""
@inline point_alias(x::AbstractFramePoint) = point_id(x)
point_alias(x::Int) = x

""" 
    @point(name, id, type=nothing)

Define a new point instance to alias the given NAIFID `id`. This macro creates an 
[`AbstractFramePoint`](@ref) subtype and its singleton instance called `name`. Its type name 
is obtained by appending `Point` to either `name` or `type` (if provided).

### Examples

```julia-repl
julia> @point Venus 299

julia> typeof(Venus)
VenusPoint 

julia> point_alias(Venus)
299

julia> @point EMB 3 EarthMoonBarycenter

julia> typeof(EMB) 
EarthMoonBarycenterPoint

julia> point_alias(EMB) 
3 
```

### See also 
See also [`@axes`](@ref) and [`point_alias`](@ref).
"""
macro point(name::Symbol, id::Int, type::Union{Symbol,Nothing}=nothing)
    # construct type name if not assigned 

    type = isnothing(type) ? name : type
    type = Symbol(format_camelcase(Symbol, String(type)), :Point)
    typ_str = String(type)
    name_str = String(name)

    pointid_expr = :(@inline Frames.point_id(::$type) = $id)
    name_expr = :(Frames.point_name(::$type) = Symbol($name_str))

    return quote
        """
            $($typ_str) <: AbstractFramePoint

        A type representing a point with NAIF ID $($id). 
        """
        struct $(esc(type)) <: AbstractFramePoint end

        """
            $($name_str)

        The singleton instance of the [`$($typ_str)`](@ref) type.
        """
        const $(esc(name)) = $(esc(type))()

        $(esc(pointid_expr))
        $(esc(name_expr))
        nothing
    end
end

""" 
    build_point(frames, name, NAIFId, class, axesid, funs; parentid, offset)

Create and add a [`FramePointNode`](@ref) to `frames` based on the input parameters. 
Current supported point classes are: `:RootPoint`, `:TimePoint`, `:EphemerisPoint`, `:FixedPoint`
and `:UpdatablePoint`.

### Inputs 
- `frames` -- Target frame system 
- `name` -- Point name, must be unique within `frames` 
- `NAIFId` -- Point NAIF ID, must be unique within `frames`
- `class` -- Point class. 
- `axesid` -- ID of the axes in which the state vector of the point is expressed. 
- `funs` -- `FramePointFunctions` object storing the functions to update the state 
            vectors of the point. It must match the type and order of `frames`

### Keywords  
- `parentid` -- NAIF ID of the parent point. Not required only for the root point.
- `offset` -- Position offset with respect to a parent point. Required only for FixedPoints.

!!! warning 
    This is a low-level function and is NOT meant to be directly used. Instead, to add a point 
    to the frame system, see [`add_point_ephemeris!`](@ref), [`add_point_fixed!`](@ref), etc...
"""
function build_point(
    frames::FrameSystem{O,T},
    name::Symbol,
    NAIFId::Int,
    class::Symbol,
    axesid::Int,
    funs::FramePointFunctions{T,O};
    parentid=nothing,
    offset=nothing,
) where {O,T}
    if has_point(frames, NAIFId)
        # Check if a point with the same NAIFId is already registered 
        # within the given FrameSystem 
        throw(
            ArgumentError(
                "A point with NAIFID $NAIFId is already registered in the input frame system.",
            ),
        )
    end

    # Check if a point with the same name does not already exist 
    if name in map(x -> x.name, frames_points(frames).nodes)
        throw(
            ArgumentError(
                "A point with name=$name is already registed in the input frame system"
            ),
        )
    end

    # Check if the given axes are known in the FrameSystem
    if !has_axes(frames, axesid)
        throw(
            ArgumentError(
                "Axes with ID $axesid are not registered in the input frame system"
            ),
        )
    end

    if isnothing(parentid)
        # If a root-point exists, check that a parent has been specified 
        if !isempty(frames_points(frames))
            throw(
                ArgumentError(
                    "A parent point is required because the input frame system " *
                    "already contains a root-point.",
                ),
            )
        end

        parentid = NAIFId # Root-point has parentid = NAIFId

    else
        # Check that the parent point is registered in frames 
        if !has_point(frames, parentid)
            throw(
                ArgumentError(
                    "The specified parent point with NAIFID $parentid is not " *
                    "registered in the input frame system.",
                ),
            )
        end
    end

    # Error check temporarily removed to avoid possible issues with unavailable point 
    # data at the programmatic start time. 

    # Check that the given functions have the correct signature 
    # for (i, fun) in enumerate((f, δf, δ²f))
    #     otype = typeof(fun(MVector{9}(zeros(T, 9)), T(1)))
    #     !(otype <: Nothing) && throw(ArgumentError(
    #         "$fun return type is $(typeof(otype)) but should be Nothing."))
    # end

    # Initialize struct caches 
    @inbounds if class in (:RootPoint, :FixedPoint)
        nzo = [@MVector zeros(Int, 2)]
        epochs = DiffCache(T[])
        stvs = [DiffCache(@MVector zeros(T, 3O))]

        if class == :FixedPoint
            for i in 1:3
                stvs[1].du[i] = offset[i]
            end
        end
    else
        # This is to handle generic frames in a multi-threading architecture 
        # without having to copy the FrameSystem
        nth = Threads.nthreads()
        nzo = [-@MVector ones(Int, 2) for _ in 1:nth]

        epochs = DiffCache(zeros(T, nth))
        stvs = [DiffCache(@MVector zeros(T, 3O)) for _ in 1:nth]
    end

    # Creates point node 
    pnode = FramePointNode{O,T,3 * O}(
        name, class, axesid, parentid, NAIFId, stvs, epochs, nzo, funs
    )

    # Insert new point in the graph
    add_point!(frames, pnode)

    # Connect the new point to the parent point in the graph 
    !isnothing(parentid) && add_edge!(frames_points(frames), parentid, NAIFId)

    return nothing
end

""" 
    add_point_root!(frames::FramesSystem, point::AbstractFramePoint, axes)

Add `point` as a root point to `frames` to initialize the points graph. Only after the 
addition of a root point, other points may be added aswell. This point is intended as the
origin, i.e., its position will equal (0., 0., 0.). `axes` is the ID or [`AbstractFrameAxes`](@ref)
instance in which the point state-vector is expressed.

!!! note
    This operation can be performed only once per [`FrameSystem`](@ref) object: multiple root 
    points in the same graph are both inadmissible and meaningless.

----

    add_point_root!(frames::FrameSystem, name::Symbol, pointid::Int, axes)

Low-level function to add point `name` with ID `pointid` to `frames` as a root-point without 
requiring the creation of an [`AbstractFramePoint`](@ref) type via the [`@point`](@ref) macro.
`axes` is the ID or [`AbstractFrameAxes`](@ref) instance in which the point state-vector is 
expressed.

### Examples 
```julia-repl
julia> FRAMES = FrameSystem{2, Float64}();

julia> @axes ICRF 1 InternationalCelestialReferenceFrame

julia> add_axes_inertial!(FRAMES, ICRF)

julia> @point SSB 0 SolarSystemBarycenter 

julia> add_point_root!(FRAMES, SSB, ICRF)

julia> @point Sun 10

julia> add_point_root!(FRAMES, Sun, ICRF)
ERROR: ArgumentError: A root-point is already registed in the given FrameSystem.
[...]
```

### See also 
See also [`add_point_ephemeris!`](@ref), [`add_point_fixed!`](@ref), [`add_point_dynamical!`](@ref)
and [`add_point_updatable!`](@ref)
"""
function add_point_root!(frames::FrameSystem, point::AbstractFramePoint, axes)
    return add_point_root!(frames, point_name(point), point_id(point), axes)
end

# Low-level function
function add_point_root!(frames::FrameSystem{O,T}, name::Symbol, pointid::Int, axes) where {O,T}

    # Check for root-point existence 
    if !isempty(frames_points(frames))
        throw(ArgumentError("A root-point is already registed in the input frame system."))
    end

    return build_point(
        frames, name, pointid, :RootPoint, axes_alias(axes), FramePointFunctions{T,O}()
    )

end


""" 
    add_point_ephemeris!(frames::FrameSystem, point::AbstractFramePoint)

Add `point` as an ephemeris point to `frames`. This function is intended for points whose 
state-vector is read from ephemeris kernels (i.e., de440.bsp). The parent point is 
automatically assigned to the point with respect to which the ephemeris data is written 
in the kernels. If that point is not available as an ephemeris point in the frame system, 
an error is thrown.

!!! note 
    The axes in which the state-vector is expressed are taken from the ephemeris data: an 
    error is returned if the axes ID is yet to be added to `frames`.

!!! warning 
    It is expected that the NAIF ID and the axes ID assigned by the user are aligned with 
    those used to generate the ephemeris kernels. No check are performed on whether these IDs
    represent the same physical bodies and axes that are intended in the kernels.

----

    add_point_ephemeris!(frames::FrameSystem, name::Symbol, naifid::Int)

Low-level function to add point `name` with ID `id` to `frames` as an ephemeris point without 
requiring the creation of an [`AbstractFramePoint`](@ref) type via the [`@point`](@ref) macro.

### Examples 
```julia-repl
julia> using Ephemerides 

julia> eph = EphemerisProvider(DE440_KERNEL_PATH);

julia> FRAMES = FrameSystem{2, Float64}(eph);

julia> @axes ICRF 1 InternationalCelestialReferenceFrame

julia> add_axes_inertial!(FRAMES, ICRF)

julia> @point SSB 0 SolarSystemBarycenter

julia> @point Sun 10 

julia> add_point_root!(FRAMES, SSB, ICRF)

julia> add_point_ephemeris!(FRAMES, Sun)

julia> @point Jupiter 599

julia> add_point_ephemeris!(FRAMES, Jupiter)
ERROR: ArgumentError: Ephemeris data for NAIFID 599 is not available in the kernels loaded in the given FrameSystem.
[...]
```

### See also 
See also [`add_point_root!`](@ref), [`add_point_fixed!`](@ref), [`add_point_dynamical!`](@ref)
and [`add_point_updatable!`](@ref)
"""
@inline function add_point_ephemeris!(frames::FrameSystem, point::AbstractFramePoint)
    return add_point_ephemeris!(frames, point_name(point), point_alias(point))
end

# Low-level function
function add_point_ephemeris!(frames::FrameSystem{O,T}, name::Symbol, naifid::Int) where {O,T}

    # Check that the kernels contain the ephemeris data for the given naifid
    if !(naifid in ephemeris_points(frames))
        throw(
            ErrorException(
                "Ephemeris data for NAIFID $naifid is not available" *
                " in the kernels loaded in the input frame system.",
            ),
        )
    end

    # Retrievet he ephemerides position records (i.e., the segment descriptors)
    pos_records = ephem_position_records(frames.eph)

    # Retrieve the parent point and the point axes from the ephemeris data 
    parentid, parent_found = 0, false
    axesid, axes_found = 0, false

    for pr in pos_records
        if pr.target == naifid
            
            # Update the parent point ID
            if !parent_found
                parentid = pr.center
                parent_found = true 

            elseif parentid != pr.center
                throw(
                    ErrorException(
                        "UnambiguityError: at least two set of data with different" * 
                        " centers are available for point with NAIFID $naifid.",
                    ),
                )

            end
            
            # Update the reference axes ID
            if !axes_found
                axesid = pr.axes 
                axes_found = true 

            elseif axesid != pr.axes 
                throw(
                    ErrorException(
                        "UnambiguityError: at least two set of data with different axes" *
                        " are available for point with NAIFID $NAIFId.",
                    ),
                )
            end

        end
    end

    # Check that the default parent is available in the FrameSystem
    if !has_point(frames, parentid)
        throw(
            ErrorException(
                "Ephemeris data for point with NAIFID $naifid is available with respect" * 
                " to point with NAIFID $parentid, which has not yet been defined" *
                " in the input frame system.",
            ),
        )
    end

    # The parent point is automatically inferred from the ephemeris kernels so it will 
    # always have available ephemeris data!

    # Checks if the axes are known to the frame system. 
    # This check is also performed by build_point, but it is reported here because 
    # it provides more specific information for ephemeris points 
    if !has_axes(frames, axesid)
        throw(
            ErrorException(
                "Ephemeris data for point with NAIFID $naifid is expressed in a set" *
                " of axes with ID $axesid, which are yet to be defined in the" *
                " input frame system.",
            ),
        )
    end

    funs = FramePointFunctions{T,O}(
        (y, t) -> ephem_compute!(y, frames.eph, DJ2000, t / DAY2SEC, naifid, parentid, 0),
        (y, t) -> ephem_compute!(y, frames.eph, DJ2000, t / DAY2SEC, naifid, parentid, 1),
        (y, t) -> ephem_compute!(y, frames.eph, DJ2000, t / DAY2SEC, naifid, parentid, 2),
        (y, t) -> ephem_compute!(y, frames.eph, DJ2000, t / DAY2SEC, naifid, parentid, 3),
    )

    return build_point(
        frames, name, naifid, :EphemerisPoint, axesid, funs; parentid=parentid
    )
end


"""
    add_point_fixed!(frames, point::AbstractFramePoint, parent, axes, offset::AbstractVector)

Add `point` as a fixed point to `frames`. Fixed points are those whose positions have a 
constant `offset` with respect their `parent` points in the given set of `axes`. Thus, points 
eligible for this class must have null velocity and acceleration with respect to `parent`.
`axes` is the ID or [`AbstractFrameAxes`](@ref) instance in which the point state-vector is 
expressed.

----

    add_point_fixed!(frames, name::Symbol, pointid::Int, parentid::Int, axes, offset::AbstractVector)

Low-level function to add point `name` with ID `id` to `frames` as a fixed-point with respect 
to point `parentid` without requiring the creation of an [`AbstractFramePoint`](@ref) type 
via the [`@point`](@ref) macro. `axes` is the ID or [`AbstractFrameAxes`](@ref) instance in 
which the point state-vector is expressed.

### Examples 
```julia-repl
julia> FRAMES = FrameSystem{2, Float64}();

julia> @axes SF -3000 SatelliteFrame

julia> add_axes_inertial!(FRAMES, SF)

julia> @point SC -10000 Spacecraft

julia> @point SolarArrayCenter -10001

julia> add_point_root!(FRAMES, SC, SF)

julia> sa_offset = [0.10, 0.15, 0.30];

julia> add_point_fixed!(FRAMES, SolarArrayCenter, SC, SF, sa_offset)
```

### See also 
See also [`add_point_root!`](@ref), [`add_point_ephemeris!`](@ref), 
[`add_point_dynamical!`](@ref) and [`add_point_updatable!`](@ref)
"""
@inline function add_point_fixed!(
    frames::FrameSystem, point::AbstractFramePoint, parent, axes, offset::AbstractVector,
)
    return add_point_fixed!(
        frames, point_name(point), point_alias(point), 
        point_alias(parent), axes_alias(axes), offset,
    )
end

# Low-level function
function add_point_fixed!(
    frames::FrameSystem{O,T}, 
    name::Symbol,
    pointid::Int, 
    parentid::Int, 
    axes, 
    offset::AbstractVector,
) where {O,T}

    if length(offset) != 3
        throw(
            DimensionMismatch(
                "The offset vector should have length 3, but has $(length(offset))."
            ),
        )
    end

    return build_point(
        frames, name, pointid, :FixedPoint, axes_alias(axes), 
        FramePointFunctions{T,O}(); parentid=parentid, offset=offset,
    )
    
end



"""
    add_point_updatable!(frames, point::AbstractFramePoint, parent, axes)

Add `point` as an updatable point to `frames`. Differently from all the other classes, the 
state vector for updatable points (expressed in the set of input `axes`) must be manually 
updated before being used for other computations. `axes` is the ID or 
[`AbstractFrameAxes`](@ref) instance in which the point state-vector is expressed. 
    
!!! note 
    This class of points becomes particularly useful if the state vector is not known a-priori, 
    e.g., when it is the output of an optimisation process which exploits the frame system.

----

    add_point_updatable!(frames, name::Symbol, pointid::Int, parentid::Int, axes)

Low-level function to add point `name` with ID `id` to `frames` as an updatable point without 
requiring the creation of an [`AbstractFramePoint`](@ref) type via the [`@point`](@ref) macro.
`axes` is the ID or [`AbstractFrameAxes`](@ref) instance in which the point state-vector is 
expressed.

### Examples 
```julia-repl
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
ERROR: UpdatablePoint with NAIFId = 1 has not been updated at time 0.1 for order 2
[...]
```

### See also 
See also [`update_point!`](@ref), [`add_point_root!`](@ref), [`add_point_ephemeris!`](@ref), 
[`add_point_dynamical!`](@ref) and [`add_point_fixed!`](@ref)
"""
function add_point_updatable!(frames::FrameSystem, point::AbstractFramePoint, parent, axes)
    return add_point_updatable!(
        frames, point_name(point), point_id(point), point_alias(parent), axes
    )
end

# Low-level function
function add_point_updatable!(
    frames::FrameSystem{O,T}, name::Symbol, pointid::Int, parentid::Int, axes
) where {O,T}

    return build_point(
        frames, name, pointid, :UpdatablePoint, axes_alias(axes), 
        FramePointFunctions{T,O}(); parentid=parentid,
    )
    
end


""" 
    add_point_dynamical!(frames, point::AbstractFramePoint, parent, axes, fun, δfun=nothing, 
        δ²fun=nothing, δ³fun=nothing)

Add `point` as a time point to `frames`. The state vector for these points depends only on 
time and is computed through the custom functions provided by the user. `axes` is the ID or 
[`AbstractFrameAxes`](@ref) instance in which the point state-vector is expressed.

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

----

    add_point_dynamical!(frames, name::Symbol, pointid::Int, parentid::Int, axes, fun, 
        δfun=nothing, δ²fun=nothing, δ³fun=nothing)

Low-level function to add point `name` with ID `id` to `frames` as a dynamical point without 
requiring the creation of an [`AbstractFramePoint`](@ref) type via the [`@point`](@ref) macro.
`axes` is the ID or [`AbstractFrameAxes`](@ref) instance in which the point state-vector is 
expressed.

### Examples 
```julia-repl
julia> FRAMES = FrameSystem{2, Float64}()

julia> @axes ICRF 1 

julia> add_axes_inertial!(FRAMES, ICRF)

julia> @point Origin 0 

julia> add_point_root!(FRAMES, Origin, ICRF)

julia> @point Satellite 1 

julia> satellite_pos(t::T) where T = [cos(t), sin(t), 0]

julia> add_point_dynamical!(FRAMES, Satellite, Origin, ICRF, satellite_pos)

julia> vector6(FRAMES, Origin, Satellite, ICRF, π/6)
6-element SVector{6, Float64} with indices SOneTo(6):
  0.8660254037844387
  0.49999999999999994
  0.0
 -0.49999999999999994
  0.8660254037844387
  0.0
```

### See also 
See also [`add_point_root!`](@ref), [`add_point_ephemeris!`](@ref),[`add_point_fixed!`](@ref)
and [`add_point_updatable!`](@ref)
"""
function add_point_dynamical!(
    frames::FrameSystem,
    point::AbstractFramePoint,
    parent,
    axes,
    fun,
    δfun = nothing,
    δ²fun = nothing,
    δ³fun = nothing,
)
    return add_point_dynamical!(
        frames, point_name(point), point_alias(point), point_alias(parent), 
        axes_alias(axes), fun, δfun, δ²fun, δ³fun,
    )
end

# Low-level function
function add_point_dynamical!(
    frames::FrameSystem{O,T},
    name::Symbol,
    pointid::Int,
    parentid::Int,
    axes,
    fun,
    δfun = nothing,
    δ²fun = nothing,
    δ³fun = nothing,
) where {O,T}

    for (order, fcn) in enumerate([δfun, δ²fun, δ³fun])
        if (O < order + 1 && !isnothing(fcn))
            @warn "ignoring $fcn, frame system order is less than $(order+1)"
        end
    end

    funs = FramePointFunctions{T,O}(
        (y, t) -> _tpoint_fun!(y, t, fun),

        # First derivative
        if isnothing(δfun)
            (y, t) -> _tpoint_δfun_ad!(y, t, fun)
        else
            (y, t) -> _tpoint_δfun!(y, t, δfun)
        end,

        # Second derivative
        if isnothing(δ²fun)
            (
                if isnothing(δfun)
                    (y, t) -> _tpoint_δ²fun_ad!(y, t, fun)
                else
                    (y, t) -> _tpoint_δ²fun_ad!(y, t, fun, δfun)
                end
            )
        else
            (y, t) -> _tpoint_δ²fun!(y, t, δ²fun)
        end,

        # Third derivative 
        if isnothing(δ³fun)
            (
                if isnothing(δ²fun)
                    (
                        if isnothing(δfun)
                            (y, t) -> _tpoint_δ³fun_ad!(y, t, fun)
                        else
                            (y, t) -> _tpoint_δ³fun_ad!(y, t, fun, δfun)
                        end
                    )
                else
                    (y, t) -> _tpoint_δ³fun_ad!(y, t, fun, δfun, δ²fun)
                end
            )
        else
            (y, t) -> _tpoint_δ³fun!(y, t, δ³fun)
        end,
    )

    return build_point(
        frames, name, pointid, :DynamicalPoint, axes_alias(axes), funs; parentid=parentid
    )

end


# Default function wrappers for time point functions! 
for (i, fun) in enumerate([:_tpoint_fun!, :_tpoint_δfun!, :_tpoint_δ²fun!, :_tpoint_δ³fun!])
    @eval begin
        function ($fun)(y, t, fn)
            @inbounds y[1:(3 * $i)] .= fn(t)
            return nothing
        end
    end
end

# Function wrapper for time-point function derivative! 
@inbounds function _tpoint_δfun_ad!(y, t, fun)
    y[1:3] .= fun(t)
    y[4:6] .= D¹(fun, t)
    return nothing
end

# Function wrappers for time-point second order derivative! 
@inbounds function _tpoint_δ²fun_ad!(y, t, fun)
    y[1:3] .= fun(t)
    y[4:6] .= D¹(fun, t)
    y[7:9] .= D²(fun, t)
    return nothing
end

@inbounds function _tpoint_δ²fun_ad!(y, t, fun, δfun)
    y[1:6] .= δfun(t)
    y[7:9] .= D²(fun, t)
    return nothing
end

# Function wrappers for time-point third order derivative! 
@inbounds function _tpoint_δ³fun_ad!(y, t, fun)
    y[1:3] .= fun(t)
    y[4:6] .= D¹(fun, t)
    y[7:9] .= D²(fun, t)
    y[10:12] .= D³(fun, t)
    return nothing
end

@inbounds function _tpoint_δ³fun_ad!(y, t, fun, δfun)
    y[1:6] .= δfun(t)
    y[7:12] .= D²(δfun, t)
    return nothing
end

@inbounds function _tpoint_δ³fun_ad!(y, t, fun, δfun, δ²fun)
    y[1:9] .= δ²fun(t)
    y[10:12] .= D³(fun, t)
    return nothing
end
