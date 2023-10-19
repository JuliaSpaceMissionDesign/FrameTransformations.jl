export add_axes_topocentric!, add_point_surface!

"""
    add_axes_topocentric!(frames, axes, λ::Number, ϕ::Number, type::Symbol, parent)

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

### See also 
See also [`add_axes_fixedoffset!`](@ref) and [`add_point_surface!`](@ref).
"""
function add_axes_topocentric!(
    frames::FrameSystem, axes::AbstractFrameAxes, λ::Number, ϕ::Number, type::Symbol, parent
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

    return add_axes_fixedoffset!(frames, axes, parent, dcm)
end

"""
    add_point_surface!(frames, point, parent, axes, pck, λ, ϕ, h::Number=0.0)

Add `point` to `frames` as a fixed point on the surface of the `parent` point. The 
relative position is specified by the longitude `λ`, the geodetic latitude `ϕ` and the altitude 
over the surface of the reference ellipsoid `h`, which defaults to 0.0. The ellipsoid 
parameters are extracted from the input TPC kernel `pck` using the NAIFId associated to the 
`parent` point.

!!! warning 
    `axes` must be a set of body-fixed reference axes for the body represented by `parent`. 
    When this constraint is not satisfied, the results may be fundamentally wrong. 

### See also 
See also [`add_point_fixed!`](@ref), [`add_axes_topocentric!`](@ref) and [`geod2pos`](@ref).
"""
function add_point_surface!(
    frames::FrameSystem,
    point::AbstractFramePoint,
    parent,
    axes,
    pck::AbstractDict,
    λ::Number,
    ϕ::Number,
    h::Number=0.0,
)
    NAIFId = point_alias(point)

    # Extract from pck constants the body parameters 
    if !haskey(pck, NAIFId) || !haskey(pck[NAIFId], :radii)
        throw(
            ArgumentError(
                "The PCK data does not contain enough information on body $NAIFId."
            ),
        )
    end

    # Compute flattening and body radius
    a, _, c = pck[NAIFId][:radii]
    f = (a - c) / a

    return add_point_surface!(frames, point, parent, axes, λ, ϕ, a, f, h)
end

"""
    add_point_surface!(frames, point, parent, axes, λ, ϕ, R, f::Number=0.0, h::Number=0.0)

Add `point` to `frames` as a fixed point on the surface of the `parent` point body. The relative 
position is specified by the longitude `λ`, the geodetic latitude `ϕ`, the reference radius 
of the ellipsoid `R` and its flattening `f`. The altitude over the reference surface of the 
ellipsoid `h` defaults to 0. 

----

add_point_surface!(frames::FrameSystem, name::Symbol, pointid::Int, parentid::Int, axesid::Int, 
    λ::Number, ϕ::Number, R::Number, f::Number=0.0, h::Number=0.0,)

    
"""
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
    pos = geod2pos(h, λ, ϕ, R, f)
    return add_point_fixed!(frames, name, pointid, parentid, axesid, pos)
end

function add_point_surface!(
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
    pos = geod2pos(h, λ, ϕ, R, f)
    return add_point_fixed!(frames, point, parent, axes, pos)
end
