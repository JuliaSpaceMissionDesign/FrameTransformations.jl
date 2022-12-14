
abstract type AbstractAxes end 
abstract type InertialAxes <: AbstractAxes end 
abstract type FixedAxes <: AbstractAxes end 
abstract type RotatingAxes <: AbstractAxes end 
abstract type ComputableAxes <: AbstractAxes end 

const _axes_classes = Dict{DataType, Int}(
    InertialAxes=>0, FixedAxes=>1, RotatingAxes=>2, 
    ComputableAxes=>3)

_get_inertial_dcm(::T, ::AbstractVector{T}) where T = DCM(T(1.0)I) 
_get_null_dcm(::T, ::AbstractVector{T}) where T = DCM(T(0.0)I)

"""
    @build_axes(frame, name, id, suptype, fun, dfun, args...)
"""
macro build_axes(frame, name::Symbol, id::Int, suptype::Symbol, 
                 fun, dfun, args...)

    frmsys = eval(frame) 
    if !(frmsys isa FrameSystem)
        throw(ArgumentError("$frame must be a FrameSystem."))
    end

    # Check that a set of axes with the same ID has not been already 
    # registered within the given FrameSystem 
    if has_vertex(axes_graph(frmsys), id) 
        throw(ErrorException("$frame already contains a set of axes with ID $id."))
    end

    # Creates a default type name 
    name_str = String(name) 
    suptyp_str = String(suptype)

    typ_str = name_str*suptyp_str 
    typ = Symbol(typ_str)

    # Only the root frame can have parent_id = id 
    parent_id = id
    parent_sym = nothing 
    parent = nothing 
    dcm = nothing 

    ref_point = 0;

    for a in args 
        a isa Expr || continue 
        a.head == :(=) || continue 

        if a.args[1] == :type 
            if !(a.args[2] isa Symbol) 
                throw(ArgumentError("Invalid type argument: $a."))
            end 

            typ = a.args[2] 
            typ_str = String(typ)
            
        elseif a.args[1] == :parent 
            if !(a.args[2] isa Symbol) 
                throw(ArgumentError("Invalid type argument: $a."))
            end 

            # Check that a set of root-axes is present 
            if isempty(axes_graph(frmsys))
                throw(ErrorException("$frame requires a set of root axes."))
            end
            
            parent_sym = a.args[2] 
            parent = eval(parent_sym)         
            parent_id = frameid(parent) 

            # Check that the parent axes belong to the given frame system 
            if !haskey(axes_graph(frmsys).ids, parent_id) 
                throw(ErrorException("The axes $(a.args[2]) are not registered in $frame."))
            end 

        elseif a.args[1] == :dcm 
            # Check correct dcm type 
            dcm = eval(a.args[2])
            if !(dcm isa DCM)
                throw(ArgumentError("Invalid type argument: $dcm. Expected DCM{T}."))
            end

        elseif a.args[1] == :refpoint 
            if !(a.args[2] isa Union{Int, Symbol})
                throw(ArgumentError("Invalid type argument: $a."))
            end

            ref_point = a.args[2] isa Int ? a.args[2] : naifid(eval(a.args[2]))
            
            # Check that the point is known within the frame system 
            if !haskey(points_graph(frmsys).ids, ref_point)
                throw(ErrorException("Point $ref_point is not registered in $frame."))
            end

            # Check that the point has as axes the same parent axes used here 
            ref_axes = get_node_from_id(points_graph(frmsys), ref_point).axes
            if ref_axes != parent_id 
                throw(ErrorException(
                    "Point $ref_point must be defined in the same axes that "*
                    "have been specified as parent for the current set. Found "*
                    "$ref_axes, expected $parent_id."))
            end

        else
            throw(ArgumentError("Unsupported keyword: $(a.args[1])."))
        end

    end

    # Check that the top frame does not already exist 
    if isnothing(parent) && !isempty(axes_graph(frmsys))
        throw(ErrorException("A set of parent axes is required because the root axes "*
                             "have already been specified."))
    end

    # Retrieves axes class
    subtyp = eval(suptype)
    class = _axes_classes[subtyp]

    # Parents\DCM checks 
    if isnothing(parent) 
        if !isnothing(dcm) # qui c'è anche la funzione da aggiungere TODO! 
            throw(ArgumentError("DCM for Root Axis is meaningless.")) 
        end 

    elseif isnothing(dcm) && class in (0, 1)
        # For non-root constant orientation axes the DCM is mandatory 
        throw(ArgumentError("The DCM from $parent_sym to $name must be specified.")) 
    end

    # Retrieves inner type of AstroAxes 
    T = get_datatype(frmsys)

    # Check that the given functions have the correct signature
    for fcn in (fun, dfun)
        otype = eval(fcn)(T(1), MVector{6}(zeros(T, 6)))
        if !(otype isa DCM{T})
            throw(ArgumentError("$fcn return type is $(typeof(otype)) "*
                                "but should be DCM{$T}."))
        end
    end

    # DCM\epoch vector settings 
    @inbounds if class in (0, 1)
        epochs = Vector{T}(undef, 0)

        R = Vector{DCM{T}}(undef, 1)
        R[1] = !isnothing(dcm) ? dcm : R[1]

        # dR is null for inertial and fixed frames
        dR = similar(R)
        dR[1] = DCM{T}(T(0)I)

    else 
        nth = Threads.nthreads()
        epochs = Vector{T}(undef, nth)
        R = Vector{DCM{T}}(undef, nth)
        dR = similar(R)
    end 

    # Checks whether the given frame is inertial 
    inertial = class == 0 || (class == 1 && isinertial(parent))

    inertial_expr = :(isinertial(::$typ) = $inertial)

    # Frame-ID Function definition 
    frameid_expr = :(frameid(::$typ) = $id)

    # Parent expression 
    parent_expr = :(parent_axes(::$typ) = $parent)

    # Create and add AstroAxes instance to AxesGraph
    reg = :(add_vertex!($(axes_graph(frmsys)), 
            AstroAxes{$T}($(QuoteNode(name)), $class, $id, $parent_id, 
                          $ref_point, $R, $dR, $epochs, $fun, $dfun)))

    # Eventually connect the axes to their parent 
    connect = !isnothing(parent) ? 
                :(connect!($(axes_graph(frmsys)), $parent_id, $id)) : :()

    doc_str = isnothing(parent) ? "root axes" : 
                "axes defined with respect to the `$(typeof(parent))`"

    return quote 
        """
            $($typ_str) <: $($suptype)
        A type representing a set of $($doc_str). 
        """
        struct $(esc(typ)) <: $(esc(suptype)) end 

        """
            $($name_str)
        The singleton instance of the [`$($typ_str)`](@ref) type.
        """
        const $(esc(name)) = $(esc(typ))()

        $(esc(frameid_expr))
        $(esc(parent_expr))
        $(esc(inertial_expr))
        $(esc(reg))
        $(esc(connect))

        nothing 
    end

end

"""
    @inertial_axes(frame, name, id, args...)
"""
macro inertial_axes(frame, name::Symbol, id::Int, args...)

    build_expr = :(@build_axes($frame, $name, $id, InertialAxes, 
                    _get_inertial_dcm, _get_null_dcm))

    for a in args 
        a isa Expr || continue 
        a.head == :(=) || continue 

        # Check for admissible parent-child relationship
        if a.args[1] == :parent 
            a.args[2] isa Symbol || throw(ArgumentError("Invalid type argument: $a"))
            
            parent = eval(a.args[2])
            if !isinertial(parent)
                # Note: if the axes have been defined as FixedAxes but are inertial, 
                # they are accepted as parents of other inertial axes! 
                throw(ErrorException("Inertial axes can only have inertial parent axes."))
            end
        end
        push!(build_expr.args, a)
    end

    return quote 
        $(esc(build_expr))
    end

end

"""
    @fixed_axes(frame, name, id, parent, dcm, args..)
"""
macro fixed_axes(frame, name::Symbol, id::Int, parent::Symbol, 
                 dcm::Union{Expr, Symbol}, args...)
    
    build_expr = :(
        @build_axes($frame, $name, $id, FixedAxes, _get_inertial_dcm,
                    _get_null_dcm, parent=$parent, dcm=$dcm))

    for a in args 
        a isa Expr || continue 
        a.head == :(=) || continue 

        push!(build_expr.args, a)
    end

    return quote 
        $(esc(build_expr))
    end
end

""" 
    @rotating_axes(frame, name, id, parent, fun, dfun, args...)
"""
macro rotating_axes(frame, name::Symbol, id::Int, parent::Symbol, 
                    fun::Symbol, dfun::Symbol, args...)

    build_expr = :(@build_axes($frame, $name, $id, RotatingAxes, $fun, 
                               $dfun, parent=$parent))
    
    for a in args 
        a isa Expr || continue 
        a.head == :(=) || continue 

        push!(build_expr.args, a)
    end

    return quote 
        $(esc(build_expr))
    end
end

"""
    @computable_axes(frame, name, id, parent::Symbol, refpoint, fun, dfun, )
"""
macro computable_axes(frame, name::Symbol, id::Int, parent::Symbol, 
                      refpoint, fun::Symbol, dfun::Symbol, args...)

    # expand in the future to allow pre-defined set of functions! 
    build_expr = :(@build_axes($frame, $name, $id, ComputableAxes, $fun, 
                            $dfun, parent=$parent, refpoint=$refpoint))

    for a in args 
        a isa Expr || continue 
        a.head == :(=) || continue 

        push!(build_expr.args, a)
    end

    return quote 
        $(esc(build_expr))
    end         

end
