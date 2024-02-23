export add_axes_topocentric!, add_point_surface!

"""
    add_axes_topocentric!(frames, axes::AbstractFrameAxes, parent, λ::Number, ϕ::Number, type::Symbol)

Add `axes` as a set of fixed-offset topocentric axes to `frames`. The orientation relative 
to the parent axes `parent` is defined throuh the longitude `λ`, the geodetic latitude `ϕ` 
and the type `type`, which may be any of the following: 

- **:NED** (North, East, Down): the X-axis points North, the Y-axis is directed eastward and 
    the Z-axis points inwards towards the nadir.
- **:SEZ** (South, East, Zenith): the X-axis points South, the Y-axis is directed East, and 
    the Z-axis points outwards towards the zenith.
- **:ENU** (East, North, Up): the X-axis points East, the Y-axis is directed North and the 
    Z-axis points outwards towards the zenith. 

!!! warning 
    The `parent` axes must be a set of body-fixed reference axes. When this constraint is 
    not satisfied, the results may be fundamentally wrong. 

----

    add_axes_topocentric!(frames, name::Symbol, axesid::Int, parentid::Int, λ, ϕ, type)

Low-level function to avoid requiring the creation of an [`AbstractFrameAxes`](@ref) via 
the [`@axes`](@ref) macro.

### See also 
See also [`add_axes_fixedoffset!`](@ref) and [`add_point_surface!`](@ref).
"""
@inline function add_axes_topocentric!(
    frames::FrameSystem, axes::AbstractFrameAxes, parent, λ::Number, ϕ::Number, type::Symbol
)
    return add_axes_topocentric!(
        frames, axes_name(axes), axes_id(axes), axes_alias(parent), λ, ϕ, type
    )
end

# Low-level function
function add_axes_topocentric!(
    frames::FrameSystem, 
    name::Symbol, 
    axesid::Int, 
    parentid::Int, 
    λ::Number,
    ϕ::Number, 
    type::Symbol
)
    if type == :NED
        dcm = angle_to_dcm(λ, -ϕ - π / 2, :ZY)
    elseif type == :SEZ
        dcm = angle_to_dcm(λ, π / 2 - ϕ, :ZY)
    elseif type == :ENU
        dcm = angle_to_dcm(λ + π / 2, π / 2 - ϕ, :ZX)
    else
        throw(ArgumentError("$type is not a supported topocentric orientation."))
    end

    return add_axes_fixedoffset!(frames, name, axesid, parentid, dcm)
end

"""
    add_point_surface!(frames, point::AbstractFramePoint, parent, axes, λ::Number, 
        ϕ::Number, R::Number, f::Number=0.0, h::Number=0.0)

Add `point` to `frames` as a fixed point on the surface of the `parent` point body. The relative 
position is specified by the longitude `λ`, the geodetic latitude `ϕ`, the reference radius 
of the ellipsoid `R` and its flattening `f`. The altitude over the reference surface of the 
ellipsoid `h` defaults to 0. 

!!! warning 
    `axes` must be a set of body-fixed reference axes for the body represented by `parent`. 
    When this constraint is not satisfied, the results may be fundamentally wrong. 

----

    add_point_surface!(frames, name::Symbol, pointid::Int, parentid::Int, axesid::Int, 
        λ::Number, ϕ::Number, R::Number, f::Number=0.0, h::Number=0.0,)

Low-level function to avoid requiring the creation of an [`AbstractFramePoint`](@ref) via 
the [`@point`](@ref) macro.
    
### See also 
See also [`add_point_fixed!`](@ref) and [`add_axes_topocentric!`](@ref).
"""
@inline function add_point_surface!(
    frames::FrameSystem,
    point::AbstractFramePoint,
    parent,
    axes,
    λ::Number,
    ϕ::Number,
    R::Number,
    f::Number=0.0,
    h::Number=0.0,
)
    return add_point_surface!(
        frames, point_name(point), point_id(point), point_alias(parent), axes_alias(axes),
        λ, ϕ, R, f, h
    )

end

@fastmath function _geod2pos(h::Number, λ::Number, ϕ::Number, R::Number, f::Number)
    # Get eccentricity from flattening 
    e² = (2 - f) * f

    sϕ, cϕ = sincos(ϕ)
    sλ, cλ = sincos(λ)

    d = R / sqrt(1 - e² * sϕ^2)
    c = (d + h) * cϕ
    s = (1 - e²) * d

    return SA[c * cλ, c * sλ, (s + h) * sϕ]
end

# Low-level function
function add_point_surface!(
    frames::FrameSystem, 
    name::Symbol,
    pointid::Int, 
    parentid::Int, 
    axesid::Int, 
    λ::Number,
    ϕ::Number,
    R::Number,
    f::Number=0.0,
    h::Number=0.0,
)
    pos = _geod2pos(h, λ, ϕ, R, f)
    return add_point_fixed!(frames, name, pointid, parentid, axesid, pos)
end

