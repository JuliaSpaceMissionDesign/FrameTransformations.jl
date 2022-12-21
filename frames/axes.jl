using Basic.Utils: format_camelcase


axes_alias(x::AbstractFrameAxes) = axes_id(x)
axes_alias(x::Int) = x

macro axes(name::Symbol, id::Int, type::Union{Symbol, Nothing}=nothing)
    # construct type name if not assigned 

    type = isnothing(type) ? name : type     
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

function build_axes(frames::FrameSystem{T}, name::Symbol, id::Int, class::Symbol, 
            f::Function, δf::Function, δ²f::Function; parentid=nothing, dcm=nothing, 
            cax_prop=ComputableAxesProperties()) where {T}

    if has_axes(frames, id)
        # Check if a set of axes with the same ID is already registered within 
        # the given frame system 
        throw(ErrorException(
            "Axes with ID = $id are already registered in the given FrameSystem."))
    end

    if name in map(x->x.name, frames_axes(frames).nodes) 
        # Check if axes with the same name also does not already exist
        throw(ErrorException(
            "Axes with name = $name are already registered in the given FrameSystem."))
    end    

    # if the frame has a parent
    if !isnothing(parentid)
        # Check if the root axes is not present
        isempty(frames_axes(frames)) && throw(ErrorException("Missing root axes."))
        
        # Check if the parent axes are registered in frame 
        if !has_axes(frames, parentid)
            throw(ErrorException("The specified parent axes with ID = $parentid are not "*
                "registered in the given FrameSystem."))
        end

    elseif class == :InertialAxes 
        parentid = id 
    end

    # Check that the given functions have the correct signature 
    for fun in (f, δf, δ²f)
        otype = fun(T(1), SVector{3}(zeros(T, 3)), SVector{3}(zeros(T, 3)))

        !(otype isa Rotation{3, T}) && throw(ArgumentError(
            "$fun return type is $(typeof(otype)) but should be Rotation{3, $T}."))
    end
    
    # Initialize struct caches
    @inbounds if class in (:InertialAxes, :FixedOffsetAxes)
        nzo = Int[]
        epochs = T[]
        R = [!isnothing(dcm) ? Rotation(dcm, DCM(T(0)I), DCM(T(0)I)) : Rotation{3}(T(1)I)]
    else
        # This is to handle generic frames in a multi-threading architecture 
        # without having to copy the FrameSystem
        nth = Threads.nthreads() 
        nzo = -ones(Int, nth)
        
        epochs = zeros(T, nth)
        R = [Rotation{3}(T(1)I) for _ = 1:nth]
    end
    
    # Creates axes node
    axnode = FrameAxesNode{T}(name, class, id, parentid, cax_prop, 
                R, epochs, nzo, f, δf, δ²f)

    # Insert the new axes in the graph
    add_axes!(frames, axnode)

    # Connect the new axes to the parent axes in the graph 
    !isnothing(parentid) && add_edge!(frames_axes(frames), parentid, id)

    nothing
end
    
_get_fixedrot9(::T, x, y) where T = Rotation{3}(T(1)I)       

function add_axes_inertial!(frames::FrameSystem{T}, axes::AbstractFrameAxes; 
            parent=nothing, dcm::Union{Nothing, DCM{T}}=nothing) where T

    name = axes_name(axes)

    # Checks for root-axes existence 
    if isnothing(parent)
        !isempty(frames_axes(frames)) && throw(ErrorException(
            "A set of parent axes for $name is required because the root axes "*
            "have already been specified in the given FrameSystem."))

        !isnothing(dcm) && throw(ArgumentError(
            "Providing a DCM for root axes is meaningless."))
    else 
        isnothing(dcm) && throw(ArgumentError(
            "Missing DCM from axes $parent."))
    end

    pid = isnothing(parent) ? nothing : axes_alias(parent)

    # construct the axes and insert in the FrameSystem
    build_axes(frames, name, axes_id(axes), :InertialAxes, 
        _get_fixedrot9, _get_fixedrot9, _get_fixedrot9; parentid=pid, dcm=dcm)

end


function add_axes_fixedoffset!(frames::FrameSystem{T}, axes::AbstractFrameAxes, 
            parent, dcm::DCM{T}) where T

    build_axes(frames, axes_name(axes), axes_id(axes), :FixedOffsetAxes, 
        _get_fixedrot9, _get_fixedrot9, _get_fixedrot9; parentid=axes_alias(parent), dcm=dcm)
end

# TODO: add computable axes 
# TODO: add rotating axes 
# TODO: add iau_axes 