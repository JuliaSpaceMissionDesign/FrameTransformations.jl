# Points

point_alias(x::AbstractFramePoint) = point_id(x)
point_alias(x::Int) = x 

macro point(name::Symbol, id::Int, type::Symbol)
    # construct type name if not assigned 

    type = Symbol(format_camelcase(Symbol, String(type)), :Point)
    typ_str = String(type)
    name_str = String(name)

    pointid_expr = :(@inline point_id(::$type) = $id)
    name_expr = :(point_name(::$type) = Symbol($name_str))

    return quote 
        """
            $($typ_str) <: AbstractFramePoint

        A type representing a point with ID $($id). 
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

function build_point(frames::FrameSystem{T, E}, name::Symbol, NAIFId::Int, class::Symbol, 
                axesid::Int, f::Function, δf::Function, δ²f::Function; 
                parentid=nothing, offset=nothing) where {T, E}

    if has_point(frames, NAIFId) 
        # Check if a point with the same NAIFId is already registered 
        # within the given FrameSystem 
        throw(ErrorException(
            "A point with NAIFId = $NAIFId is already registered in the given FrameSystem."))
    end

    # Check if a point with the same name does not already exist 
    if name in map(x->x.name, frames_points(frames).nodes)
        throw(ErrorException(
            "A point with name = $name is already registed in the given FrameSystem"))
    end 

    # Check if the given axes are known in the FrameSystem
    !has_axes(frames, axesid) && throw(ErrorException(
            "Axes with ID = $axesid are not registered in the given FrameSystem"))
            
    if isnothing(parentid) 
        # If a root-point exists, check that a parent has been specified 
        if !isempty(frames_points(frames)) 
            throw(ErrorException("A parent point is required because the given FrameSystem "*
                "already contains a root-point."))
        end

        parentid = NAIFId # Root-point has parentid = NAIFId

    else 
        # Check that the parent point is registered in frames 
        if !has_point(frames, parentid)
            throw(ErrorException("The specified parent point with NAIFId = $parentid is not "*
                "registered in the given FrameSystem"))
        end
    end

    # Check that the given functions have the correct signature 
    for (i, fun) in enumerate((f, δf, δ²f))
        # TO DO: da cambiare con l'epooca minima per le effemeridi! 
        otype = typeof(fun(MVector{9}(zeros(T, 9)), T(1)))
        !(otype <: Nothing) && throw(ArgumentError(
            "$fun return type is $(typeof(otype)) but should be Nothing."))
    end

    # Initialize struct caches 
    @inbounds if class in (:RootPoint, :FixedPoint)
        nzo = Int[]
        epochs = T[]
        stvs = [@MVector zeros(T, 9)]

        if class == :FixedPoint
            for i = 1:3 
                stvs[i] = offset[i]
            end
        end
    else 
        # This is to handle generic frames in a multi-threading architecture 
        # without having to copy the FrameSystem
        nth = Threads.nthreads() 
        nzo = -ones(Int, nth)
        
        epochs = zeros(T, 9)
        stvs = [@MVector zeros(T, 9) for _ = 1:nth]
    end

    # Creates point node 
    pnode = FramePointNode{T}(name, class, axesid, parentid, NAIFId, 
                stvs, epochs, nzo, f, δf, δ²f)

    # Insert new point in the graph
    add_point!(frames, pnode)

    # Connect the new point to the parent point in the graph 
    !isnothing(parentid) && add_edge!(frames_points(frames), parentid, NAIFId)

    nothing 
end


_empty_stv_update!(::AbstractVector{T}, ::T) where {T} = nothing


function add_point_root!(frames::FrameSystem, name::Symbol, NAIFId::Int, axes)

    # Check for root-point existence 
    if !isempty(frames_points(frames)) 
        throw(ErrorException(
            "A root-point is already registed in the given FrameSystem."))
    end

    build_point(frames, name, NAIFId, :RootPoint, axes_alias(axes), 
                _empty_stv_update!, _empty_stv_update!, _empty_stv_update!)

end


function add_point_ephemeris!(frames::FrameSystem, name::Symbol, NAIFId::Int, parent)

    # Check that the kernels contain the ephemeris data for the given NAIFId
    if !(NAIFId in ephemeris_points(frames))
        throw(ErrorException("Ephemeris data for NAIFId = $NAIFId is not available "*
            "in the kernels loaded in the given FrameSystem."))
    end

    # Check that the parent point is admissible
    parentid = point_alias(parent) 
    parentclass = get_node(frames_points(frames), parentid).class
    if !(parentclass in (:RootPoint, :EphemerisPoint))
        throw(ErrorException("The specified parent point with NAIFId = $parentid is a "*
            "$parentclass in the given FrameSystem, but only RootPoints and "*
            "EphemerisPoints are accepted as parents of EphemerisPoints."))
    end

    # Check that the parent point has available ephemeris data 
    if !(parentid in ephemeris_points(frames)) 
        throw(ErrorException("Insufficient ephemeris data has been loaded to compute "*
            "the point with NAIFId = $NAIFId with respect to the parent point with "*
            "NAIFId = $parentid"))
    end

    # Retrieves the axes stored in the ephemeris kernels for the given point
    axesid = nothing 
    for pr in ephem_position_records(frames.eph)
        if pr.target == NAIFId
            if isnothing(axesid)
                axesid = pr.frame 
            elseif axesid != pr.frame 
                throw(ErrorException("UnambiguityError: at least two set of data "*
                    "with different axes are available for point with NAIFId = $NAIFId."))
            end
        end
    end

    # Checks if the axes are known to the frame system. 
    # This check is also performed by build_point, but it is reported here because 
    # it provides more specific information for ephemeris points 
    if !has_axes(frames, axesid)
        throw(ErrorException("Ephemeris data for point with NAIFId = $NAIFId is expressed "*
            "in a set of axes with ID = $axesid, which are yet to be defined in the "*
            "given FrameSystem."))
    end

    build_point(frames, name, NAIFId, :EphemerisPoint, axesid, 
                (y, t) -> ephem_compute_order!(y, frames.eph, t, 0., NAIFId, parentid, 0),
                (y, t) -> ephem_compute_order!(y, frames.eph, t, 0., NAIFId, parentid, 1),
                (y, t) -> ephem_compute_order!(y, frames.eph, t, 0., NAIFId, parentid, 2),; 
                parentid=parentid)

end 


function add_point_fixed!(frames::FrameSystem{T, E}, name::Symbol, NAIFId::Int, parent, 
            axes, offset::AbstractVector{T}) where {T, E}

    
    if length(offset) != 3
        throw(DimensionMismatch(
            "The offset vector has length 3, but has $(length(offset))."))
    end

    build_point(frames, name, NAIFId, :FixedPoint, axes_alias(axes), 
                _empty_stv_update!, _empty_stv_update!, _empty_stv_update!; 
                parentid=point_alias(parent), offset=offset)
        
end


function add_point_updatable!(frames::FrameSystem{T, E}, name::Symbol, NAIFId::Int, 
                              parent, axes) where {T, E}

    build_point(frames, name, NAIFId, :UpdatablePoint, axes_alias(axes), 
                _empty_stv_update!, _empty_stv_update!, _empty_stv_update!; 
                parentid=point_alias(parent))
end


@inline function add_point_root!(frames::FrameSystem, point::AbstractFramePoint, axes) 
    add_point_root!(frames, point_name(point), point_id(point), axes)    
end

@inline function add_point_ephemeris!(frames::FrameSystem, point::AbstractFramePoint, parent)
    add_point_ephemeris!(frames, point_name(point), point_id(point), parent)
end

@inline function add_point_fixed!(frames::FrameSystem, point::AbstractFramePoint, 
                    parent, axes, offset) 
    add_point_fixed!(frames, point_name(point), point_id(point), parent, axes, offset)    
end

@inline function add_point_updatable!(frames, point::AbstractFramePoint, parent, axes) 
    add_point_updatable!(frames, point_name(point), point_id(point), parent, axes)    
end

