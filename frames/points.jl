
abstract type AbstractPoint end 
abstract type RootPoint <: AbstractPoint end  
abstract type FixedPoint <: AbstractPoint end 
abstract type TimePoint <: AbstractPoint end 
abstract type EphemerisPoint <: TimePoint end  
abstract type UpdatablePoint <: AbstractPoint end # da modificare a mano! 

const _point_classes = Dict{DataType, Int}(
    RootPoint=>0, TimePoint=>1, EphemerisPoint=>2, 
    FixedPoint=>3, UpdatablePoint=>4)

_empty_stv_update!(::AbstractVector{T}, ::T) where {T} = nothing

@inline function update_point!(frame::FrameSystem, point::AbstractPoint, 
                                y::AbstractVector)
    update_point!(frame, naifid(point), y)
end

# Updates updatable points 
function update_point!(frame::FrameSystem{T}, naifid::Int, 
                        y::AbstractVector{T}, ep::Number) where T
    
    point = get_node_from_id(points_graph(frame), naifid)
    
    tid = Threads.threadid()
    @inbounds if point.epochs[tid] != ep 
        point.epochs[tid] = ep
        @simd for i in eachindex(y)
            point.stv[tid][i] = y[i]
        end
    end
end


"""
    @build_point(frame, name, id, suptype, fun, axes, args...)
"""
macro build_point(frame, name::Symbol, id::Int, subtype::Symbol,
                  fun, dfun, axes::Union{Int, Symbol}, args...)

    # missing check on fun signature!

    frmsys = eval(frame)
    if !(frmsys isa FrameSystem)
        throw(ArgumentError("$frame must be a FrameSystem object."))
    end 

    # Check that a point with the same NAIFId has not been already 
    # registerd within the given FrameSystem 
    if has_vertex(points_graph(frmsys), id)
        throw(ErrorException("$frame already contains a point with NAIFId $id."))
    end 

    # Check that the given axes are known in the frame system
    axes_id = axes isa Int ? axes : frameid(eval(axes))
    if !haskey(axes_graph(frmsys).ids, axes_id)
        throw(ErrorException("Axes $axes are not registered in $frame."))
    end

    # Creates a default type name 
    name_str = String(name)
    subtype_str = String(subtype)

    typ_str = name_str*subtype_str; 
    typ = Symbol(typ_str)

    parent_id = 0 
    parent_sym = nothing 
    parent = nothing 
    
    offset = nothing

    # Retrieves inner type of AstroPoint 
    T = get_datatype(frmsys);
    
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

            # Check that a root-point is present 
            if isempty(points_graph(frmsys))
                throw(ErrorException("$frmsys requires a root point."))
            end

            parent_sym = a.args[2]
            parent = eval(parent_sym)
            parent_id = naifid(parent)

            # Check that the parent point belongs to the given frame system 
            if !haskey(points_graph(frmsys).ids, parent_id)
                throw(ErrorException("The point $(a.args[2]) is not registered in $frame."))
            end
        elseif a.args[1] == :offset
            offset = eval(a.args[2])
            if !(offset isa AbstractVector && length(offset) == 3)
                throw(ArgumentError("Invalid type argument: $offset. Expected a Vector{$T} of length 3."))
            end 

        else 
            throw(ArgumentError("Unsupported keyword: $(a.args[1])"))
        end
    end

    # If a root-point exists, check that a parent has been specified!
    if isnothing(parent) && !isempty(points_graph(frmsys))
        throw(ErrorException("A parent point is required!"))
    end 

    # Point class
    class = _point_classes[eval(subtype)]

    # State vector\epoch settings 
    if class == 0 || class == 3

        epochs = Vector{T}(undef, 0)
        stvs = Vector{MVector{6, T}}(undef, 1)

        # ocio a sto zeros che potrebbe dare problemi con tipi strani
        @inbounds stvs[1] = @MVector zeros(T, 6)

        if class == 3 
            @inbounds for i = 1:3
                stvs[1][i] = offset[i]
            end
        end
    else 
        nthreads = Threads.nthreads() 

        epochs = Vector{T}(undef, nthreads)
        stvs = Vector{MVector{6, T}}(undef, nthreads)

        for i in eachindex(stvs)
            stvs[i] = @MVector zeros(T, 6)
        end

    end

    # Check if ephemeris data is available for that point 
    hasephe = id in available_ephemeris_bodies(frmsys)

    # NAIF-ID function definition 
    naifid_expr = :(naifid(::$typ) = $id)

    # Parent expression
    parent_expr = :(parent_point(::$typ) = $parent)

    # # Create and add AstroPoint instance to PointSystem 
    reg = :(add_vertex!($(points_graph(frmsys)), 
            AstroPoint{$T}($(QuoteNode(name)), $hasephe, $class, $axes_id,
                           $parent_id, $id, $stvs, $epochs, $fun, $dfun)))

    # Eventually connect the point to its parent 
    connect = !isnothing(parent) ? 
                :(connect!($(points_graph(frmsys)), $parent_id, $id)) : :()

    doc_str = isnothing(parent) ? "root point" :
                " point defined with respect to `$(typeof(parent))`"
    return quote 
        """ 
            $($typ_str) <: $($subtype)
        A type representing a $($doc_str).
        """
        struct $(esc(typ)) <: $(esc(subtype)) end 

        """
            $($name_str)
        The singleton instance of the [`$($typ_str)`](@ref) type.  
        """
        const $(esc(name)) = $(esc(typ))()

        $(esc(naifid_expr))
        $(esc(parent_expr))

        # Registers the point 
        $(esc(reg))
        $(esc(connect))

        nothing
    end

end


"""
    @root_point(frame, name, naifid, axes, args...)
"""
macro root_point(frame, name::Symbol, naifid::Int, 
                 axes::Union{Int, Symbol}, args...)

    build_expr = :(@build_point($frame, $name, $naifid, RootPoint, 
                            _empty_stv_update!, _empty_stv_update!, $axes))

    for a in args 
        a isa Expr || continue 
        a.head == :(=) || continue 

        if a.args[1] == :parent 
            throw(ArgumentError("Root points cannot have parent points!"))        
        end 

        push!(build_expr.args, a)
    end

    return quote 
        $(esc(build_expr))
    end

end


"""
    @ephemeris_point(frame, name, naifid, parent, args...)
"""
macro ephemeris_point(frame, name::Symbol, NAIFId::Int, 
                      parent::Symbol, args...)

    frmsys = eval(frame)
    if !(frmsys isa FrameSystem)
        throw(ArgumentError("$frame must be a FrameSystem object."))
    end

    # Check that the kernels contain the data for naifid 
    if !(NAIFId in available_ephemeris_bodies(frmsys))
        throw(ErrorException("Ephemeris data for point $NAIFId is not available"*
                             " in the kernels loaded in $(frame)!"))
    end

    # Check that the parent point is admissible and has available ephemeris data 
    parent_id = naifid(eval(parent))
    parent_point = points_graph(frmsys).nodes[points_graph(frmsys).ids[parent_id]]
    if !parent_point.heph
        throw(ErrorException("Insufficient ephemeris data has been loaded to compute "*
                             "the point $NAIFId with respect to point $(parent_id)!"))

    elseif !(parent_point.class in (0, 2))
        throw(ErrorException("Only RootPoints and other EphemerisPoints are "*
                            "accepted as parents of EphemerisPoints."))
    end

    # Retrieves the axes
    axes_id = nothing
    for pr in get_position_records(frmsys.eph)
        if pr.target == NAIFId 
            if isnothing(axes_id)
                axes_id = pr.frame
            elseif axes_id != pr.frame 
                throw(ErrorException("Unambiguity error: at least two set of data"*
                      " with different axes are available for point $(NAIFId)!"))
            end
        end
    end

    # This check is also performed by @build_point but is reported 
    # here because it provide more specific information when the axes 
    # are unknown within the frame system
    if !haskey(axes_graph(frmsys).ids, axes_id)
        throw(ErrorException("Ephemeris data for point $NAIFId is expressed in a "*
              "set of axes with ID $axes_id, which are yet to be defined in $(frame)!"))
    end

    fun = (y, t) -> unsafe_compute_order!(y, frmsys.eph, t, 0., NAIFId, parent_id, 
                                            cph_use_naifid+cph_sec+cph_km, 0)

    dfun = (y, t) -> unsafe_compute_order!(y, frmsys.eph, t, 0., NAIFId, parent_id, 
                                            cph_use_naifid+cph_sec+cph_km, 1)

    build_expr = :(@build_point($frame, $name, $NAIFId, EphemerisPoint, 
                                $fun, $dfun, $axes_id, parent=$parent))
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
    @fixed_point(frame, name, naifid, parent, axes, offset, args...)
"""
macro fixed_point(frame, name::Symbol, NAIFId::Int, parent::Symbol, 
                  axes::Union{Int, Symbol}, 
                  offset::Union{Expr, Symbol}, args...)

    build_expr = :(@build_point($frame, $name, $NAIFId, FixedPoint, 
                    _empty_stv_update!, _empty_stv_update!, 
                    $axes, offset=$offset, parent=$parent))

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
    @time_point(frame, name, naifid, parent, axes, fun, dfun, args...)
"""
macro time_point(frame, name::Symbol, NAIFId::Int, parent::Symbol,
                 axes::Union{Int, Symbol}, fun::Symbol, dfun::Symbol, args...)

    build_expr = :(@build_point($frame, $name, $NAIFId, TimePoint, 
                    $fun, $dfun, $axes, parent=$parent))

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
    @updatable_point(frame, name, naifid, parent, axes, args...)
"""
macro updatable_point(frame, name::Symbol, NAIFId::Int, parent::Symbol,
                      axes::Union{Int, Symbol}, args...)

    build_expr = :(@build_point($frame, $name, $NAIFId, UpdatablePoint, 
                    _empty_stv_update!, _empty_stv_update!, 
                    $axes, parent=$parent))

    for a in args 
        a isa Expr || continue 
        a.head == :(=) || continue 

        push!(build_expr.args, a)
    end

    return quote 
        $(esc(build_expr))
    end
end

