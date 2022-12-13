import ForwardDiff.derivative 

abstract type AbstractAxes end 
abstract type InertialAxes <: AbstractAxes end 
abstract type FixedAxes <: AbstractAxes end 
abstract type RotatingAxes <: AbstractAxes end 
abstract type ComputableAxes <: AbstractAxes end 


_get_fixrot_1(::T, x, y) where T = Rotation{T}(DCM(T(1)I)) 
_get_fixrot_2(::T, x, y) where T = Rotation{T}(DCM(T(1)I), DCM(T(0)I)) 
_get_fixrot_3(::T, x, y) where T = Rotation{T}(DCM(T(1)I), DCM(T(0)I), DCM(T(0)I)) 


"""
    @build_axes(frame, name, id, class, fun, dfun, ddfun, args...)
"""
macro build_axes(frm, name::Symbol, id::Int, class, 
                 fun, dfun, ddfun, args...)

    frame = eval(frm) 
    if !(frame isa FrameSystem)
        throw(ArgumentError("$frm must be a FrameSystem."))
    end

    # Retrieves inner type of FrameAxesNode  
    T = get_datatype(frame)

    # Check that a set of axes with the same ID has not been already 
    # registered within the given FrameSystem 
    if has_axes(frame, id)
        throw(ErrorException("$frm already contains a set of axes with ID $id."))
    end

    # Creates a default type name 
    name_str = String(name) 
    suptyp_str = String(class)

    typ_str = name_str*suptyp_str 
    typ = Symbol(typ_str)

    parent_id = id    # Only the root frame can have parent_id = id 
    parent = nothing 
    dcm = nothing 

    # ref_point = 0;

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
            
            parent = eval(a.args[2])         
            parent_id = axesid(parent) 

            # Check that a set of root-axes is present 
            if isempty(axes_graph(frame))
                throw(ErrorException("$frm requires a set of root axes."))
            end

            # Check that the parent axes belong to the given frame system 
            if !has_axes(frame, parent_id)
                throw(ErrorException("The axes $(a.args[2]) are not registered in $frm."))
            end 

        elseif a.args[1] == :dcm 

            # Check correct dcm type 
            dcm = eval(a.args[2])
            if !(dcm isa DCM)
                throw(ArgumentError("Invalid type argument: $dcm. Expected DCM{$T}."))
            end

        # elseif a.args[1] == :refpoint 
        #     if !(a.args[2] isa Union{Int, Symbol})
        #         throw(ArgumentError("Invalid type argument: $a."))
        #     end

        #     ref_point = a.args[2] isa Int ? a.args[2] : naifid(eval(a.args[2]))
            
        #     # Check that the point is known within the frame system 
        #     if !haskey(points_graph(frmsys).ids, ref_point)
        #         throw(ErrorException("Point $ref_point is not registered in $frame."))
        #     end

        #     # Check that the point has as axes the same parent axes used here 
        #     ref_axes = get_node_from_id(points_graph(frmsys), ref_point).axes
        #     if ref_axes != parent_id 
        #         throw(ErrorException(
        #             "Point $ref_point must be defined in the same axes that "*
        #             "have been specified as parent for the current set. Found "*
        #             "$ref_axes, expected $parent_id."))
        #     end

        else
            throw(ArgumentError("Unsupported keyword: $(a.args[1])."))
        end

    end

    # Check that the top frame does not already exist 
    if isnothing(parent) && !isempty(axes_graph(frame))
        throw(ErrorException("A set of parent axes is required because the root axes "*
                             "have already been specified."))
    end

    # Parents\DCM checks 
    if isnothing(parent) 
        if !isnothing(dcm) # qui c'è anche la funzione da aggiungere TODO! 
            throw(ArgumentError("DCM for Root Axis is meaningless.")) 
        end 

    elseif isnothing(dcm) && class in (:InertialAxes, :FixedAxes)
        # For non-root constant orientation axes the DCM is mandatory 
        throw(ArgumentError("The DCM from $parent to $name must be specified.")) 
    end

    # Check that the given functions have the correct signature
    for (i, fcn) in enumerate([fun, dfun, ddfun])
        otype = eval(fcn)(T(1), SVector{3}(zeros(T, 3)), 
                                SVector{3}(zeros(T, 3)))

        if !(otype isa Rotation{i, T})
            throw(ArgumentError("$fcn return type is $(typeof(otype)) "*
                                "but should be Rotation{$i, $T}."))
        end
    end

    # Rotation\epoch vector settings 
    @inbounds if class in (:InertialAxes, :FixedAxes)
        epochs = Vector{T}(undef, 0)

        R = Vector{Rotation{3, T}}(undef, 1)
        R[1] = !isnothing(dcm) ? Rotation{T}(dcm, DCM(T(0)I), DCM(T(0)I)) : R[1]
        nzo = -ones(Int, 1)

        inertial = class == :InertialAxes || isinertial(parent)

    else 
        nth = Threads.nthreads()
        epochs = Vector{T}(undef, nth)
        R = Vector{Rotation{3, T}}(undef, nth)
        nzo = -ones(Int, nth)

        inertial = false 
    end 

    # Sets default computable axes 
    if class != :ComputableAxes 
        comp_axes_prop = ComputableAxesProperties()
    end

    # Is inertial expression
    inertial_expr = :(isinertial(::$typ) = $inertial)

    # Axes-ID Function definition 
    axesid_expr = :(axesid(::$typ) = $id)

    # Parent expression 
    parent_expr = :(parent_axes(::$typ) = $parent)

    # Create and add AstroAxes instance to AxesGraph
    reg = :(add_vertex!($(axes_graph(frame)), 
            FrameAxesNode{$T}($(QuoteNode(name)), Symbol($class), $id, 
                        $parent_id, $comp_axes_prop, $R, $epochs, $nzo, 
                          $fun, $dfun, $ddfun)))

    # Eventually connect the axes to their parent 
    connect = !isnothing(parent) ? 
                :(connect!($(axes_graph(frame)), $parent_id, $id)) : :()

    doc_str = isnothing(parent) ? "root axes" : 
                "axes defined with respect to the `$(typeof(parent))`"

    return quote 
        """
            $($typ_str) <: $($class)
        A type representing a set of $($doc_str). 
        """
        struct $(esc(typ)) <: $(esc(class)) end 

        """
            $($name_str)
        The singleton instance of the [`$($typ_str)`](@ref) type.
        """
        const $(esc(name)) = $(esc(typ))()

        $(esc(axesid_expr))
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
                    _get_fixrot_1, _get_fixrot_2, _get_fixrot_3))

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
        @build_axes($frame, $name, $id, FixedAxes, _get_fixrot_1, 
                _get_fixrot_2, _get_fixrot_3, parent=$parent, dcm=$dcm))

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
    @rotating_axes(frame, name, id, parent, fun, args...)
"""
macro rotating_axes(frame, name::Symbol, id::Int, parent::Symbol, 
                    fn::Symbol, dfun=nothing, ddfun=nothing, args...)

    fun = eval(fn)

    if isnothing(dfun) 
        dfun = t -> Rotation(fun(t), derivative(fun, t)) 
    end 

    if isnothing(ddfun) 
        ddfun = t -> Rotation(fun(t), df(t), derivative(fun, t))
    end 

    f = (t, x, y) -> fun(t)
    df = (t, x, y) -> dfun(t)
    ddf = (t, x, y) -> ddfun(t)

    build_expr = :(@build_axes($frame, $name, $id, RotatingAxes, $f, 
                               $df, $ddf, parent=$parent))
    
    for a in args 
        a isa Expr || continue 
        a.head == :(=) || continue 

        ## controlla fun e dfun 

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

