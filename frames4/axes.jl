using Basic.Utils: format_camelcase

const AXES_CLASSES = (
    :InertialAxes,
    :FixedOffsetAxes
)

axes_alias(x::AbstractFrameAxes) = axes_id(x)
axes_alias(x::Int) = x

macro axes(name::Symbol, id::Int, type::Symbol)
    # construct type name if not assigned 

    type = Symbol(format_camelcase(Symbol, String(type)), :Axes)
    typ_str = String(type)
    name_str = String(name)

    axesid_expr = :(@inline axes_id(::$type) = $id)
    name_expr = :(axes_name(::$type) = Symbol($name_str))
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

function build_axes(
    frames::FrameSystem{T, E}, name::Symbol, id::Int, class::Symbol, 
    f::Function, δf::Function, δ²f::Function; 
    parentid=nothing, dcm=nothing, cax_prop=ComputableAxesProperties()
) where {T, E}

    if has_axes(frames, id)
        # Check if a set of axes with the same ID is already registered within 
        # the given frame system 
        throw(ErrorException(
            "Axes with ID = $id are already registered in the FrameSystem.")
        )
    end

    if name in map(x->x.name, frames_axes(frames).nodes) 
        # Check if axes with the same name also does not already exists
        throw(ErrorException(
        "Axes with name=$name are already registered in the FrameSystem.")
    )
    end    

    # if the frame has a parent
    if !isnothing(parentid)
        # Check if the root axes is not present
        isempty(frames_axes(frames)) && throw(ErrorException("Missing root axes."))
        
        # Check if the parent axes are not registered in frame 
        if !has_axes(frames, parentid)
            throw(ErrorException("Parent axes with ID=$id are not registered in the FrameSystem."))
        end
    elseif class == :InertialAxes 
        # TODO: why?
        parentid = id 
    end

    # Check that the given functions have the correct signature 
    for (i, fun) in enumerate((f, δf, δ²f))
        otype = fun(T(1), SVector{3}(zeros(T, 3)), SVector{3}(zeros(T, 3)))

        !(otype isa Rotation{i, T}) && throw(ArgumentError(
            "$fun return type is $(typeof(otype)) but should be Rotation{$i, $T}."))
    end

    # Initialize struct caches
    @inbounds if class in (:InertialAxes, :FixedOffsetAxes)
        nzo = -ones(Int, 1)
        epochs = Vector{T}(undef, 0)
        R = Vector{Rotation{3, T}}(undef, 1)
        R[1] = !isnothing(dcm) ? Rotation(dcm, DCM(T(0)I), DCM(T(0)I)) : R[1]
    else
        # this is to handle generic frames in a multi-threading architecture without 
        # having to copy the FrameSystem
        nth = Threads.nthreads() 
        nzo = -ones(Int, nth)
        epochs = Vector{T}(undef, nth)
        R = Vector{Rotation{3, T}}(undef, nth)
    end
    
    # Creates axes node
    axnode = FrameAxesNode{T}(
        name, class, id, parentid, cax_prop, 
        R, epochs, nzo,
        f, δf, δ²f
    )

    # Insert the new axes in the graph
    add_axes!(frames, axnode)

    # Connect the new axes to the parent axes in the graph 
    if !isnothing(parentid)
        add_edge!(frames_axes(frames), parentid, id)
    end

    nothing
end

_get_fixedrot3(::T, x, y) where T = Rotation(DCM(T(1)I)) 
_get_fixedrot6(::T, x, y) where T = Rotation(DCM(T(1)I), DCM(T(0)I)) 
_get_fixedrot9(::T, x, y) where T = Rotation(DCM(T(1)I), DCM(T(0)I), DCM(T(0)I)) 

function axes_register_inertial!(frames::FrameSystem, name::Symbol, id::Int; 
    parent=nothing, dcm=nothing)

    # Checks for root-axes existence 
    if isnothing(parent)
        !isempty(frames_axes(frames)) && throw(ErrorException(
            "A set of parent axes for $name is required because the root axes "*
            "have already been specified."))

        !isnothing(dcm) && throw(ArgumentError("Providing a DCM for root axes is meaningless."))
    else 
        isnothing(dcm) && throw(ArgumentError("Missing DCM from axes $parent."))
    end

    pid = isnothing(parent) ? nothing : axes_alias(parent)

    # construct the axes and insert in the FrameSystem
    build_axes(frames, name, id, :InertialAxes, 
        _get_fixedrot3, _get_fixedrot6, _get_fixedrot9; parentid=pid, dcm=dcm)

end

function axes_register_inertial!(frames, axes::A; args...) where {A<:AbstractFrameAxes}
    axes_register_inertial!(frames, axes_name(axes), axes_id(axes); args...)
end

function axes_register_fixedoffset!(frames::FrameSystem, name::Symbol, id::Int, parent, 
    dcm::DCM{T}) where {T}
    # construct the axes 
    build_axes(frames, name, id, :FixedOffsetAxes, 
        _get_fixedrot3, _get_fixedrot6, _get_fixedrot9; parentid=axes_alias(parent), dcm=dcm)
end

function axes_register_fixedoffset!(frames, axes::A, parent, dcm) where {A<:AbstractFrameAxes}
    axes_register_fixedoffset!(frames, axes_name(axes), axes_id(axes), parent, dcm)
end