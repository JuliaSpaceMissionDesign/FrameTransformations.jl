import ForwardDiff.derivative 

abstract type AbstractAxes end 

get_alias(x::AbstractAxes) = axesid(x)
get_alias(x::Int) = x 

_get_fixrot_1(::T, x, y) where T = Rotation(DCM(T(1)I)) 
_get_fixrot_2(::T, x, y) where T = Rotation(DCM(T(1)I), DCM(T(0)I)) 
_get_fixrot_3(::T, x, y) where T = Rotation(DCM(T(1)I), DCM(T(0)I), DCM(T(0)I)) 

macro axes(name::Symbol, id::Int, type=nothing)

    if isnothing(type)
        type = Symbol(:Abstract, name, :Axes)
    end        

    typ_str = String(type)
    name_str = String(name)

    axesid_expr = :(axesid(::$type) = $id)
    name_expr = :(axes_name(::$type) = Symbol($name_str))

    return quote 
        """
            $($typ_str) <: AbstractAxes
        A type representing a set of axes with ID $($id). 
        """
        struct $(esc(type)) <: AbstractAxes end

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


function build_axes(frame::FrameSystem{T}, 
                    name::Symbol, id::Int, class::Symbol, 
                    fun::Function, dfun::Function, ddfun::Function; 
                    parentid=nothing, # oppure Int  
                    dcm=nothing, 
                    caprop=ComputableAxesProperties()) where {T}
    
    # Check that a set of axes with the same ID is not already registered 
    # within the given frame system 
    has_axes(frame, id) && throw(ErrorException(
            "Axes with ID=$id are already registered in the given frame system."))
    
    # Check that a similar name also does not already exists
    name in map(x->x.name, axes_graph(frame).nodes) && throw(ErrorException(
        "Axes with name=$name are already registered in the given frame system."))

    if !isnothing(parentid)
        # Check that a set of root axes is present 
        isempty(axes_graph(frame)) && throw(ErrorException(
            "Missing root axes."))

        # Check that the parent axes are registered in frame 
        !has_axes(frame, parentid) && throw(ErrorException(
            "The parent axes $id are not registered in this frame system"))

    elseif class == :InertialAxes 
        parentid = id 
    end


    # Check that the given functions have the correct signature 
    for (i, fcn) in enumerate([fun, dfun, ddfun])
        otype = fcn(T(1), SVector{3}(zeros(T, 3)), SVector{3}(zeros(T, 3)))

        !(otype isa Rotation{i, T}) && throw(ArgumentError(
            "$fcn return type is $(typeof(otype)) but should be Rotation{$i, $T}."))
    end 

    # Initializations of epoch\rotation vectors 
    @inbounds if class in (:InertialAxes, :FixedAxes) 
        nzo = -ones(Int, 1)
        epochs = Vector{T}(undef, 0)

        R = Vector{Rotation{3, T}}(undef, 1)
        R[1] = !isnothing(dcm) ? Rotation(dcm, DCM(T(0)I), DCM(T(0)I)) : R[1]

    else 
        nth = Threads.nthreads() 
        nzo = -ones(Int, nth)
        epochs = Vector{T}(undef, nth)
        R = Vector{Rotation{3, T}}(undef, nth)
    end

    # IT = Tuple{T, SVector{3, T}, SVector{3, T}}

    # Creates axes object
    axesnode = FrameAxesNode{T}(name, class, id, parentid, caprop, R, 
                         epochs, nzo, fun, dfun, ddfun)

    # Adds the new axes inside the graph 
    add_axes!(frame, axesnode)

    # Connects axes to their parent axes in the graph 
    !isnothing(parentid) && connect!(axes_graph(frame), parentid, id) 

end


"""
    add_inertial_axes!(frame, id, parent, dcm, name)
"""
function add_inertial_axes!(frame::FrameSystem, name::Symbol, id::Int;
                            parent=nothing, dcm=nothing)

    # Checks for root-axes existence 
    if isnothing(parent)
        
        !isempty(axes_graph(frame)) && throw(ErrorException(
            "A set of parent axes is required because the root axes "*
            "have already been specified."))

        !isnothing(dcm) && throw(ArgumentError(
            "Providing a DCM for root axes is meaningless."))

    else 
        isnothing(dcm) && throw(ArgumentError(
            "Missing rotation DCM from axes $parent."))
    
    end

    pid = isnothing(parent) ? nothing : get_alias(parent)

    build_axes(frame, name, id, :InertialAxes, _get_fixrot_1, 
        _get_fixrot_2, _get_fixrot_3; parentid=pid, dcm=dcm)

end 


"""
    add_fixed_axes!(frame, name, id, parent, dcm)
"""
function add_fixed_axes!(frame::FrameSystem{T}, name::Symbol, id::Int, 
                         parent, dcm::DCM{T}) where T 

    build_axes(frame, name, id, :FixedAxes, _get_fixrot_1, 
            _get_fixrot_2, _get_fixrot_3; dcm=dcm, 
            parentid=get_alias(parent))

end


"""
    add_rotating_axes!(frame, name, id, parent, fun[, dfun[, ddfun]])
"""
function add_rotating_axes!(frame::FrameSystem{T}, name::Symbol, id::Int, 
                            parent, fun, dfun=nothing, ddfun=nothing) where T

    build_axes(frame, name, id, :RotatingAxes, 
                (t, x, y) -> Rotation(fun(t)), 
                isnothing(dfun) ? 
                    (t, x, y) -> Rotation(fun(t), derivative(fun, t)) : 
                    (t, x, y) -> Rotation(dfun(t)),

                isnothing(ddfun) ?
                    (isnothing(dfun) ? 
                        (t, x, y) -> Rotation(fun(t), derivative(fun, t), derivative(τ->derivative(fun, τ), t)) : 
                        (t, x, y) -> Rotation(dfun(t)..., derivative(τ->derivative(fun, τ), t))) : 
                    (t, x, y) -> Rotation(ddfun),

                parentid=get_alias(parent))

end

"""
    add_computable_axes!(frame, name, id, parent, ???)
"""
function add_computable_axes!(frame::FrameSystem, name::Symbol, id::Int, 
                              parent)

    
    

end

add_rotating_axes!(frame, axes::AbstractAxes, parent, fun, args...) = 
    add_rotating_axes(frame, axes_name(axes), axesid(axes), parent, fun, args...)

add_fixed_axes!(frame, axes::AbstractAxes, parent, dcm) = 
    add_fixed_axes!(frame, axes_name(axes), axesid(axes), parent, dcm)
    
add_inertial_axes!(frame, axes::AbstractAxes; args...) =
    add_inertial_axes!(frame, axes_name(axes), axesid(axes); args...)
