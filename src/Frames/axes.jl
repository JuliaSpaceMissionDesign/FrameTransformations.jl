export @axes,
    add_axes_inertial!,
    add_axes_rotating!,
    add_axes_fixedoffset!,
    add_axes_computable!,
    add_axes_projected!,
    add_axes_ephemeris!,
    is_inertial,
    is_timefixed,
    axes_alias

"""
    axes_name(axes::AbstractFrameAxes)

Return the name of `axes`.
"""
function axes_name end

""" 
    axes_id(axes::AbstractFrameAxes)

Return the ID associated to `axes`.
"""
function axes_id end

""" 
    is_inertial(frame::FrameSystem, axes::AbstractFrameAxes)
    is_inertial(frame::FrameSystem, axesid::Int)

Return true if the given axes are inertial, i.e., non rotating with respect to the root inertial 
axes.

!!! note 
    FixedOffsetAxes with respect to an inertial set of axes, are also consired inertial.
"""
function is_inertial(frame::FrameSystem, axes::AbstractFrameAxes)
    return is_inertial(frame, axes_alias(axes))
end
is_inertial(frame::FrameSystem, axesid::Int) = is_inertial(frames_axes(frame), axesid)

function is_inertial(axframe::MappedNodeGraph, axesid::Int)
    node = get_node(axframe, axesid)
    if node.class in (:InertialAxes, :FixedOffsetAxes, :ProjectedAxes)
        if node.id != node.parentid
            return is_inertial(axframe, node.parentid)
        else
            return true # Root axes are always inertial
        end
    else
        return false
    end
end

""" 
    is_timefixed(frame::FrameSystem, axes::AbstractFrameAxes)
    is_timefixed(frame::FrameSystem, axesid::Int)

Return true if the given axes are time-fixed, i.e., their orientation does not change in 
time with respect to the root inertial axes. 

!!! note 
    Only `:InertialAxes` and `:FixedOffsetAxes` defined with respect to other inertial axes 
    are here considered as time fixed. 
"""
function is_timefixed(frame::FrameSystem, axes::AbstractFrameAxes)
    return is_timefixed(frame, axes_alias(axes))
end
is_timefixed(frame::FrameSystem, axesid::Int) = is_timefixed(frames_axes(frame), axesid)

function is_timefixed(axframe::MappedNodeGraph, axesid::Int)
    node = get_node(axframe, axesid)
    if node.class in (:InertialAxes, :FixedOffsetAxes)
        if node.id != node.parentid
            return is_timefixed(axframe, node.parentid)
        else
            return true # Root axes are always time fixed
        end
    else
        return false
    end
end

"""
    axes_alias(ax::AbstractFrameAxes)

Return the axes ID. 
"""
axes_alias(x::AbstractFrameAxes) = axes_id(x)
axes_alias(x::Int) = x

"""
    @axes(name, id, type=nothing)

Define a new axes instance to alias the given `id`. This macro creates an 
[`AbstractFrameAxes`](@ref) subtype and its singleton instance callen `name`. Its type name 
is obtained by appending `Axes` to either `name` or `type` (if provided). 

### Examples 

```julia-repl
julia> @axes ICRF 1 InternationalCelestialReferenceFrame

julia> typeof(ICRF)
InternationalCelestialReferenceFrameAxes

julia> axes_alias(ICRF) 
1

julia> @axes IAU_EARTH 10013

julia> typeof(IAU_EARTH)
IauEarthAxes
```

### See also 
See also [`@point`](@ref) and [`axes_alias`](@ref).
"""
macro axes(name::Symbol, id::Int, type::Union{Symbol,Nothing}=nothing)
    # construct type name if not assigned 

    type = isnothing(type) ? name : type
    type = Symbol(format_camelcase(Symbol, String(type)), :Axes)
    typ_str = String(type)
    name_str = String(name)

    axesid_expr = :(@inline Frames.axes_id(::$type) = $id)
    name_expr = :(Frames.axes_name(::$type) = Symbol($name_str))

    return quote
        """
            $($typ_str) <: AbstractFrameAxes

        A type representing a set of axes with ID $($id). 
        """
        struct $(esc(type)) <: AbstractFrameAxes end

        """
            $($name_str)

        The singleton instance of the [`$($typ_str)`](@ref) type.
        """
        const $(esc(name)) = $(esc(type))()

        $(esc(axesid_expr))
        $(esc(name_expr))
        nothing
    end
end

"""
    build_axes(frames, name, id, class, funs; parentid, dcm, cax_prop)

Create and add a [`FrameAxesNode`](@ref) to `frames` based on the input parameters. Current 
supported classes are: `:InertialAxes`, `:FixedOffsetAxes`, `:RotatingAxes`, `:ProjectedAxes`, 
`:EphemerisAxes` and `:ComputableAxes`.

### Inputs 
- `frames` -- Target frame system 
- `name` -- Axes name, must be unique within `frames` 
- `id` -- Axes ID, must be unique within `frames`
- `class` -- Axes class.  
- `funs` -- `FrameAxesFunctions` object storing the functions to compute the DCM and, 
            eventually, its time derivatives. It must match the type and order of `frames`.

### Keywords 
- `parentid` -- Axes ID of the parent axes. Not required only for the root axes.
- `dcm` -- DCM with respect to the parent axes. Required only for FixedOffsetAxes. 
- `cax_prop` -- `ComputableAxesProperties`, required only by ComputableAxes.

!!! warning 
    This is a low-level function and is NOT meant to be directly used. Instead, to add a set of
    axes to the frame system, see [`add_axes_inertial!`](@ref), [`add_axes_rotating!`](@ref), 
    [`add_axes_computable!`](@ref), [`add_axes_fixedoffset!`](@ref) and [`add_axes_projected!`](@ref).

"""
function build_axes(
    frames::FrameSystem{O,T},
    name::Symbol,
    id::Int,
    class::Symbol,
    funs::FrameAxesFunctions{T,O};
    parentid=nothing,
    dcm=nothing,
    cax_prop=ComputableAxesProperties(),
) where {O,T}
    if has_axes(frames, id)
        # Check if a set of axes with the same ID is already registered within 
        # the given frame system 
        throw(
            ArgumentError(
                "Axes with ID $id are already registered in the input frame system."
            ),
        )
    end

    if name in map(x -> x.name, frames_axes(frames).nodes)
        # Check if axes with the same name also do not already exist
        throw(
            ArgumentError(
                "Axes with name=$name are already registered in the input frame system."
            ),
        )
    end

    # if the frame has a parent
    if !isnothing(parentid)
        # Check if the root axes is not present
        isempty(frames_axes(frames)) && throw(ArgumentError("Missing root axes."))

        # Check if the parent axes are registered in frame 
        if !has_axes(frames, parentid)
            throw(
                ArgumentError(
                    "The specified parent axes with ID $parentid are not " *
                    "registered in the input frame system.",
                ),
            )
        end

    elseif class == :InertialAxes
        parentid = id
    end

    # Check that the given functions have the correct signature 
    for i in 1:O
        otype = funs[i](T(1), SVector{3O}(zeros(T, 3O)), SVector{3O}(zeros(T, 3O)))

        !(otype isa Rotation{O,T}) && throw(
            ArgumentError(
                "$(funs[i]) return type is $(typeof(otype)) but should be Rotation{$O, $T}.",
            ),
        )
    end

    # Initialize struct caches
    @inbounds if class in (:InertialAxes, :FixedOffsetAxes)
        nzo = Int[]
        epochs = T[]
        angles = [DiffCache(@MVector zeros(T, 3O))]
        R = [!isnothing(dcm) ? Rotation{O}(dcm) : Rotation{O}(T(1)I)]
    
    else
        # This is to handle generic frames in a multi-threading architecture 
        # without having to copy the FrameSystem
        nth = Threads.nthreads()
        nzo = -ones(Int, nth)
        epochs = zeros(T, nth)
        R = [Rotation{O}(T(1)I) for _ in 1:nth]
        angles = [DiffCache(@MVector zeros(T, 3O)) for _ in 1:nth]

    end

    # Creates axes node
    axnode = FrameAxesNode{O,T,3 * O}(
        name, class, id, parentid, cax_prop, R, epochs, nzo, funs, angles
    )

    # Insert the new axes in the graph
    add_axes!(frames, axnode)

    # Connect the new axes to the parent axes in the graph 
    !isnothing(parentid) && add_edge!(frames_axes(frames), parentid, id)

    return nothing
end

"""
    add_axes_inertial!(frames, axes; parent=nothing, dcm=nothing)

Add `axes` as a set of inertial axes to `frames`. Only inertial axes can be used as root axes 
to initialise the axes graph. Only after the addition of a set of inertial axes, other axes 
classes may be added aswell. Once a set of root-axes has been added, `parent` and `dcm` 
become mandatory fields.

!!! note
    The parent of a set of inertial axes must also be inertial. FixedOffset axes that are 
    that only have inertial parents are also accepted.

----

    add_axes_inertial!(frames, name::Symbol, axesid::Int; parentid=nothing, dcm=nothing)

Low-level function to add axes `name` with id `axesid` to `frames` without requiring the 
creation of an [`AbstractFrameAxes`](@ref) type via the [`@axes`](@ref) macro.

### Examples 
```julia-repl 
julia> FRAMES = FrameSystem{2, Float64}();

julia> @axes ICRF 1 InternationalCelestialReferenceFrame 

julia> add_axes_inertial!(FRAMES, ICRF)

julia> @axes ECLIPJ2000 17 

julia> add_axes_inertial!(FRAMES, ECLIPJ2000)
ERROR: A set of parent axes for ECLIPJ2000 is required
[...]

julia> add_axes_inertial!(FRAMES, ECLIPJ2000; parent=ICRF, dcm=angle_to_dcm(π/3, :Z))
```

### See also 
See also [`add_axes_rotating!`](@ref), [`add_axes_fixedoffset!`](@ref) and [`add_axes_computable!`](@ref) 
"""
function add_axes_inertial!(
    frames::FrameSystem,
    axes::AbstractFrameAxes;
    parent = nothing,
    dcm = nothing,
)
    pid = isnothing(parent) ? nothing : axes_alias(parent)
    return add_axes_inertial!(frames, axes_name(axes), axes_id(axes); parentid=pid, dcm=dcm)
end

function add_axes_inertial!(
    frames::FrameSystem{O,T},
    name::Symbol,
    axesid::Int;
    parentid = nothing,
    dcm::Union{Nothing,DCM{T}} = nothing,
) where {O,T}

    # Checks for root-axes existence 
    if isnothing(parentid)
        if !isempty(frames_axes(frames))
            throw(
                ArgumentError(
                    "A set of parent axes for $name is required because the root axes " *
                    "have already been specified in the given FrameSystem.",
                ),
            )
        end

        if !isnothing(dcm)
            throw(ArgumentError("Providing a DCM for root axes is meaningless."))
        end

    else
        # Check that the parent axes are contained in frames! 
        if !has_axes(frames, parentid)
            throw(
                ArgumentError(
                    "The parent axes with ID $parentid are not present in the given FrameSystem."
                ),
            )
        end

        # Check that the parent axes are inertial\timefixed
        if !is_timefixed(frames, parentid)
            throw(ArgumentError("The parent axes for inertial axes must also be inertial."))
        end

        if isnothing(dcm)
            throw(ArgumentError("Missing DCM from axes with ID $parentid."))
        end
    end

    return build_axes(
        frames,
        name,
        axesid,
        :InertialAxes,
        FrameAxesFunctions{T,O}();
        parentid=parentid,
        dcm=dcm,
    )
end

"""
    add_axes_fixedoffset!(frames::FrameSystem, axes, parent, dcm::DCM)

Add `axes` as a set of fixed offset axes to `frames`. Fixed offset axes have a constant 
orientation with respect to their `parent` axes, represented by `dcm`, a 
Direction Cosine Matrix (DCM).

!!! note 
    While inertial axes do not rotate with respect to the star background, fixed offset axes 
    are only constant with respect to their parent axes, but might be rotating with respect 
    to some other inertial axes.

----

    add_axes_fixedoffset!(frames, name::Symbol, axesid::Int, parentid::Int, dcm:DCM)
   
Low-level function to add axes `name` with id `axesid` to `frames` with a fixed-offset from 
`parentid` without requiring the creation of an [`AbstractFrameAxes`](@ref) type via the 
[`@axes`](@ref) macro.

### Examples 
```julia-repl 
julia> FRAMES = FrameSystem{1, Float64}();

julia> @axes ICRF 1 InternationalCelestialReferenceFrame 

julia> add_axes_inertial!(FRAMES, ICRF)

julia> @axes ECLIPJ2000 17 

julia> add_axes_fixedoffset!(FRAMES, ECLIPJ2000, ICRF, angle_to_dcm(π/3, :Z))
```

### See also 
See also [`add_axes_rotating!`](@ref), [`add_axes_inertial!`](@ref) and [`add_axes_computable!`](@ref) 
"""
@inline function add_axes_fixedoffset!(
    frames::FrameSystem, axes::AbstractFrameAxes, parent, dcm::DCM
)
    return add_axes_fixedoffset!(
        frames, axes_name(axes), axes_id(axes), axes_alias(parent), dcm
    )
end

# Low-level function
function add_axes_fixedoffset!(
    frames::FrameSystem{O,T}, name::Symbol, axesid::Int, parentid::Int, dcm::DCM{T}
) where {O,T}

    return build_axes(
        frames, name, axesid, :FixedOffsetAxes, FrameAxesFunctions{T,O}(); 
        parentid = parentid, dcm = dcm
    )
end

"""
    add_axes_rotating!(frames, axes, parent, fun, δfun=nothing, δ²fun=nothing, δ³fun=nothing)

Add `axes` as a set of rotating axes to `frames`. The orientation of these axes depends only 
on time and is computed through the custom functions provided by the user. 

The input functions must accept only time as argument and their outputs must be as follows: 

- **fun**: return a Direction Cosine Matrix (DCM).
- **δfun**: return the DCM and its 1st order time derivative.
- **δ²fun**: return the DCM and its 1st and 2nd order time derivatives.
- **δ³fun**: return the DCM and its 1st, 2nd and 3rd order time derivatives.

If `δfun`, `δ²fun` or `δ³fun` are not provided, they are computed via automatic differentiation.

!!! warning 
    It is expected that the input functions and their outputs have the correct signature. This 
    function does not perform any checks on the output types. 

----

    add_axes_rotating!(frames, name::Symbol, axesid::Int, parentid::Int, fun, δfun=nothing, 
        δ²fun=nothing, δ³fun=nothing)
   
Low-level function to add axes `name` with id `axesid` to `frames` as rotating axes without 
requiring the creation of an [`AbstractFrameAxes`](@ref) type via the [`@axes`](@ref) macro.

### Examples 
```julia-repl 
julia> FRAMES = FrameSystem{3, Float64}();

julia> @axes Inertial 1

julia> add_axes_inertial!(FRAMES, Inertial)

julia> @axes Synodic 2 

julia> fun(t) = angle_to_dcm(t, :Z);

julia> add_axes_rotating!(FRAMES, Synodic, Inertial, fun)

julia> R = rotation6(FRAMES, Inertial, Synodic, π/6);

julia> R[1]
DCM{Float64}:
  0.866025  0.5       0.0
 -0.5       0.866025  0.0
  0.0       0.0       1.0

julia> R[2]
DCM{Float64}:
 -0.5        0.866025  0.0
 -0.866025  -0.5       0.0
  0.0        0.0       0.0
```

### See also 
See also [`add_axes_fixedoffset!`](@ref), [`add_axes_inertial!`](@ref) and [`add_axes_computable!`](@ref) 
"""
@inline function add_axes_rotating!(
    frames::FrameSystem,
    axes::AbstractFrameAxes,
    parent,
    fun,
    δfun = nothing,
    δ²fun = nothing,
    δ³fun = nothing
)
    return add_axes_rotating!(
        frames, axes_name(axes), axes_id(axes), axes_alias(parent), fun, δfun, δ²fun, δ³fun
    )
end

# Low-level function
function add_axes_rotating!(
    frames::FrameSystem{O,T},
    name::Symbol,
    axesid::Int,
    parentid::Int,
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

    funs = FrameAxesFunctions{T,O}(
        (t, x, y) -> Rotation{O}(fun(t)),

        # First derivative 
        if isnothing(δfun)
            (t, x, y) -> Rotation{O}(fun(t), D¹(fun, t))
        else
            (t, x, y) -> Rotation{O}(δfun(t))
        end,

        # Second derivative 
        if isnothing(δ²fun)
            (
                if isnothing(δfun)
                    (t, x, y) -> Rotation{O}(fun(t), D¹(fun, t), D²(fun, t))
                else
                    (t, x, y) -> Rotation{O}(δfun(t)..., D²(fun, t))
                end
            )
        else
            (t, x, y) -> Rotation{O}(δ²fun(t))
        end,

        # Third derivative 
        if isnothing(δ³fun)
            (
                if isnothing(δ²fun)
                    (
                        if isnothing(δfun)
                            (t, x, y) ->
                                Rotation{O}(fun(t), D¹(fun, t), D²(fun, t), D³(fun, t))
                        else
                            (t, x, y) -> Rotation{O}(δfun(t)..., D²(δfun, t)...)
                        end
                    )
                else
                    (t, x, y) -> Rotation{O}(δ²fun(t)..., D³(fun, t))
                end
            )
        else
            (t, x, y) -> Rotation{O}(δ³fun(t))
        end,
    )

    return build_axes(frames, name, axesid, :RotatingAxes, funs; parentid=parentid)
end

"""
    add_axes_computable!(frames, axes, parent, v1, v2, seq::Symbol)

Add `axes` as a set of computable axes to `frames`. Computable axes differ from rotating axes 
because they are computed through two vectors that are defined within the frame system itself. 
Computable axes are the equivalent of SPICE's parameterized two-vector frames. 

These axes are built such that the first vector, known as the primary vector, is parallel to 
one axis of the frame; the component of the secondary vector orthogonal to the first is parallel
to another axis of the frame, and the cross product of the two vectors is parallel to the 
remaining axis. 

The primary and secondary vectors, `v1` and `v2` are instances of `ComputableAxesVector`, 
which is used to define the NAIF IDs of the vector origin and target, and its order. Current 
accepted order values are: 1 (position), 2 (velocity) and 3 (acceleration). 

For example, to define a primary vector that is parallel to the Sun's (NAIF ID = 10) velocity 
with respect to the Solary system barycenter (NAIF ID = 0), `v1` must be set as: 
`v1 = ComputableAxesVector(10, 0, 1)`.

`seq` is a combination of two letters that is used to identify the desired pointing 
directions of the primary and secondary vectors. Accepted sequences are: `:XY`, `:YX`, `:XZ`,
`:ZX`, `:YZ` and `:ZY`. 

Given a spacecraft registered as a point in the frame system, an example of a set of computable 
axes is the Local Vertical Local Horizon (LVLH), where the spacecraf's nadir direction and 
velocity direction define the axes orientation.  

!!! note 
    Regardless of the original set of axes in which the primary and secondary vectors are 
    defined, the axes orientation is automatically computed by rotating them to `parent`.

----

    add_axes_computable!(frames, name::Symbol, axesid::Int, parentid::Int, v1, v2, seq)
   
Low-level function to add axes `name` with id `axesid` to `frames` as computable axes without 
requiring the creation of an [`AbstractFrameAxes`](@ref) type via the [`@axes`](@ref) macro.

### Examples 
```julia-repl 
julia> using Ephemerides 

julia> eph = EphemerisProvider(DE440_KERNEL_PATH);

julia> FRAMES = FrameSystem{4, Float64}(eph);

julia> @point SSB 0 SolarySystemBarycenter 

julia> @point Sun 10 

julia> @axes ICRF 1 InternationalCelestialReferenceFrame 

julia> add_axes_inertial!(FRAMES, ICRF)

julia> add_point_root!(FRAMES, SSB, ICRF)

julia> add_point_ephemeris!(FRAMES, Sun)

julia> @axes SunFrame 2

julia> v1 = ComputableAxesVector(10, 0, 1)
ComputableAxesVector(10, 0, 1)

julia> v2 = ComputableAxesVector(10, 0, 2)
ComputableAxesVector(10, 0, 2)

julia> add_axes_computable!(FRAMES, SunFrame, ICRF, v1, v2, :XY)
```

### See also 
See also [`ComputableAxesVector`](@ref), [`add_axes_fixedoffset!`](@ref), [`add_axes_inertial!`](@ref) 
and [`add_axes_computable!`](@ref)
"""
@inline function add_axes_computable!(
    frames::FrameSystem,
    axes::AbstractFrameAxes,
    parent,
    v1::ComputableAxesVector,
    v2::ComputableAxesVector,
    seq::Symbol,
)
    return add_axes_computable!(
        frames, axes_name(axes), axes_id(axes), axes_alias(parent), v1, v2, seq
    )
end

# Low-level function
function add_axes_computable!(
    frames::FrameSystem{O,T},
    name::Symbol,
    axesid::Int,
    parentid::Int,
    v1::ComputableAxesVector,
    v2::ComputableAxesVector,
    seq::Symbol,
) where {O,T}
    if !(seq in (:XY, :YX, :XZ, :ZX, :YZ, :ZY))
        throw(
            ArgumentError("$seq is not a valid rotation sequence for two vectors frames.")
        )
    end

    for v in (v1, v2)
        for id in (v.from, v.to)
            if !has_point(frames, id)
                throw(
                    ArgumentError(
                        "Point with NAIFID $id is unknown in the given frame system."
                    ),
                )
            end
        end
    end

    funs = FrameAxesFunctions{T,O}(
        (t, x, y) -> Rotation{O}(twovectors_to_dcm(x, y, seq)),
        (t, x, y) -> Rotation{O}(_two_vectors_to_rot6(x, y, seq)),
        (t, x, y) -> Rotation{O}(_two_vectors_to_rot9(x, y, seq)),
        (t, x, y) -> Rotation{O}(_two_vectors_to_rot12(x, y, seq)),
    )

    return build_axes(
        frames, name, axesid, :ComputableAxes, funs; 
        parentid = parentid, cax_prop = ComputableAxesProperties(v1, v2),
    )

end


"""
    add_axes_projected!(frames, axes, parent, fun)

Add `axes` as a set of projected axes to `frames`. The orientation of these axes depends only 
on time and is computed through the custom functions provided by the user. 

Projected axes are similar to rotating axis, except that all the positions, velocity, etc ... 
are rotated by the 0-order rotation (i.e. the derivatives of the rotation matrix are null, 
despite the rotation depends on time).

!!! warning 
    It is expected that the input function and their outputs have the correct signature. This 
    function does not perform any checks on the output types. 

----

    add_axes_projected!(frames, name::Symbol, axesid::Int, parentid::Int, fun)
   
Low-level function to add axes `name` with id `axesid` to `frames` as projected axes without 
requiring the creation of an [`AbstractFrameAxes`](@ref) type via the [`@axes`](@ref) macro.
"""
function add_axes_projected!(frames::FrameSystem, axes::AbstractFrameAxes, parent, fun)
    return add_axes_projected!(
        frames, axes_name(axes), axes_id(axes), axes_alias(parent), fun
    )
end

# Low-level function
function add_axes_projected!(
    frames::FrameSystem{O,T}, name::Symbol, axesid::Int, parentid::Int, fun
) where {O,T}
    funs = FrameAxesFunctions{T,O}(
        (t, x, y) -> Rotation{O}(fun(t)),
        (t, x, y) -> Rotation{O}(fun(t), DCM(0.0I)),
        (t, x, y) -> Rotation{O}(fun(t), DCM(0.0I), DCM(0.0I)),
        (t, x, y) -> Rotation{O}(fun(t), DCM(0.0I), DCM(0.0I), DCM(0.0I)),
    )

    return build_axes(frames, name, axesid, :ProjectedAxes, funs; parentid=parentid)
end



"""
    add_axes_ephemeris!(frames, axes, rot_seq::Symbol)
    add_axes_ephemeris!(frames, axes, fun, δfun=nothing, δ²fun=nothing, δ³fun=nothing)

Add `axes` as a set of ephemeris axes to `frames`. The orientation of these axes is computed 
with a series of 3 rotations specified by `rot_seq`. The euler angles and their derivatives 
needed by the rotation are extracted from the ephemeris kernels loaded in `frames`.
The parent axes are automatically assigned to the axes with respect to which the orientation 
data has been written in the kernels.

This operation is only possible if the ephemeris kernels loaded within `frames` contain 
orientation data for the AXES ID associated to `axes`. An error is returned if the parent 
axes ID is yet to be added to `frames.`

The rotation sequence is defined by a `Symbol` specifing the rotation axes. This function
assigns `dcm = A3 * A2 * A1` in which `Ai` is the DCM related with the **i-th** rotation. 
The possible `rot_seq` values are: `:XYX`, `XYZ`, `:XZX`, `:XZY`, `:YXY`, `:YXZ`, `:YZX`, 
`:YZY`, `:ZXY`, `:ZXZ`, `:ZYX`, or `:ZYZ`.

Alternatively, one can provide custom functions that take a vector of euler angles and their 
derivatives and return the DCM and its derivatives. 

The input functions must accept only time as argument and their outputs must be as follows: 

- **fun**: return a Direction Cosine Matrix (DCM).
- **δfun**: return the DCM and its 1st order time derivative.
- **δ²fun**: return the DCM and its 1st and 2nd order time derivatives.
- **δ³fun**: return the DCM and its 1st, 2nd and 3rd order time derivatives.

!!! warning 
    It is expected that the axes ID assigned by the user are aligned with those used to 
    generate the ephemeris kernels. No check are performed on whether these IDs represent 
    the same physical axes that are intended in the kernels.

----

    add_axes_ephemeris!(frames, name::Symbol, axesid::Int, rot_seq::Symbol)
    add_axes_ephemeris!(frames, name::Symbol, axesid::Int, fun, δfun=nothing, δ²fun=nothing, δ³fun=nothing)

Low-level functions to add axes `name` with id `axesid` to `frames` as ephemeris axes without 
requiring the creation of an [`AbstractFrameAxes`](@ref) type via the [`@axes`](@ref) macro.

### See also 
See also [`ComputableAxesVector`](@ref), [`add_axes_fixedoffset!`](@ref), [`add_axes_inertial!`](@ref) 
and [`add_axes_computable!`](@ref)
"""
function add_axes_ephemeris!(frames::FrameSystem, axes::AbstractFrameAxes, rot_seq::Symbol)
    return add_axes_ephemeris!(frames, axes_name(axes), axes_id(axes), rot_seq)
end

function add_axes_ephemeris!(
    frames::FrameSystem, 
    axes::AbstractFrameAxes, 
    fun, 
    δfun=nothing, 
    δ²fun=nothing, 
    δ³fun=nothing,
)
    return add_axes_ephemeris!(
        frames, axes_name(axes), axes_id(axes), fun, δfun, δ²fun, δ³fun
    )
end

# Low-level function
function add_axes_ephemeris!(
    frames::FrameSystem{O,T}, name::Symbol, axesid::Int, rot_seq::Symbol
) where {O,T}

    # Check and retrieve the parent ID for the given axes
    parentid = _check_axes_ephemeris(frames, axesid)

    if rot_seq in (:ZYX, :XYX, :XYZ, :XZX, :XZY, :YXY, :YXZ, :YZX, :YZY, :ZXY, :ZXZ, :ZYZ)
        funs = FrameAxesFunctions{T,O}(
            (t, x, y) -> Rotation{O}(_3angles_to_rot3(x, rot_seq)),
            (t, x, y) -> Rotation{O}(_3angles_to_rot6(x, rot_seq)),
            (t, x, y) -> Rotation{O}(_3angles_to_rot9(x, rot_seq)),
            (t, x, y) -> Rotation{O}(_3angles_to_rot12(x, rot_seq)),
        )
    else
        throw(ArgumentError("The rotation sequence :$rot_seq is not valid."))
    end

    return build_axes(frames, name, axesid, :EphemerisAxes, funs; parentid=parentid)

end

# Low-level function
function add_axes_ephemeris!(
    frames::FrameSystem{O,T}, 
    name::Symbol, 
    axesid::Int, 
    fun, 
    δfun = nothing, 
    δ²fun = nothing, 
    δ³fun = nothing
) where {O,T}

    # Check and retrieve the parent ID for the given axes
    parentid = _check_axes_ephemeris(frames, axesid)

    if (
        (O ≥ 2 && isnothing(δfun)) ||
        (O ≥ 3 && isnothing(δ²fun)) ||
        (O ≥ 4 && isnothing(δ³fun))
    )
        throw(
            ArgumentError(
                "Transformation requires function derivatives up to the $(O-1) order."
            ),
        )
    end

    funs = FrameAxesFunctions{T,O}(
        (t, x, y) -> Rotation{O}(fun(x)),
        (t, x, y) -> Rotation{O}(δfun(x)),
        (t, x, y) -> Rotation{O}(δ²fun(x)),
        (t, x, y) -> Rotation{O}(δ³fun(x)),
    )

    return build_axes(frames, name, axesid, :EphemerisAxes, funs; parentid=parentid)
    
end

# Perform checks on the requested ephemeris axes
function _check_axes_ephemeris(frames::FrameSystem, axesid::Int)

    # Check that the kernels contain orientation data for the given axes ID 
    if !(axesid in ephemeris_axes(frames))
        throw(
            ArgumentError(
                "Orientation data for AXESID $axesid is not available " *
                "in the kernels loaded in the given FrameSystem.",
            ),
        )
    end

    # Retrieve the parent from the ephemeris data 
    orient_records = ephem_orient_records(frames.eph)
    parentid = nothing

    for or in orient_records
        if or.target == axesid
            if isnothing(parentid)
                parentid = or.axes
            elseif parentid != or.axes
                throw(
                    ErrorException(
                        "UnambiguityError: at least two set of orientation data " *
                        "with different centers are available for axes with AXESID $NAIFId.",
                    ),
                )
            end
        end
    end

    # Check that the default parent is available in the FrameSystem! 
    if !has_axes(frames, parentid)
        throw(
            ArgumentError(
                "Orientation data for axes with AXESID $axesid is available " *
                "with respect to axes with AXESID $parentid, which has not yet been defined " *
                "in the given FrameSystem.",
            ),
        )
    end

    return parentid
end

# Functions for an easier definition\handling of ephemeris axes
@inline function _3angles_to_rot3(θ, seq::Symbol)
    @inbounds angle_to_dcm(θ[1], θ[2], θ[3], seq)
end

@inline function _3angles_to_rot6(θ, seq::Symbol)
    return _3angles_to_rot3(θ, seq), _3angles_to_δdcm(θ, seq)
end

@inline function _3angles_to_rot9(θ, seq::Symbol)
    return (_3angles_to_rot3(θ, seq), _3angles_to_δdcm(θ, seq), _3angles_to_δ²dcm(θ, seq))
end

@inline function _3angles_to_rot12(θ, seq::Symbol)
    return (
        _3angles_to_rot3(θ, seq),
        _3angles_to_δdcm(θ, seq),
        _3angles_to_δ²dcm(θ, seq),
        _3angles_to_δ³dcm(θ, seq),
    )
end
