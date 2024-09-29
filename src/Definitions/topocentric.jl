
"""
    add_axes_topocentric!(frames, name::Symbol, id::Int, parent, λ, ϕ, mount)

Add topocentric axes to `frames` at a specified location and mounting.

The orientation relative to the parent axes `parent` is defined throuh the longitude `λ`, 
the geodetic latitude `ϕ` and the mounting type `mount`, which may be any of the following: 

- `:NED` (North, East, Down): the X-axis points North, the Y-axis is directed eastward and 
    the Z-axis points inwards towards the nadir.
- `:SEZ` (South, East, Zenith): the X-axis points South, the Y-axis is directed East, and 
    the Z-axis points outwards towards the zenith.
- `:ENU` (East, North, Up): the X-axis points East, the Y-axis is directed North and the 
    Z-axis points outwards towards the zenith. 

!!! warning 
    The parent axes must be a set of body-fixed reference axes. This is under user 
    resposibility. 
"""
function add_axes_topocentric!(
    frames::FrameSystem, name::Symbol, id::Int, parent,
    λ::Number, ϕ::Number, mount::Symbol
)
    if mount == :NED
        dcm = angle_to_dcm(λ, -ϕ - π / 2, :ZY)
    elseif mount == :SEZ
        dcm = angle_to_dcm(λ, π / 2 - ϕ, :ZY)
    elseif mount == :ENU
        dcm = angle_to_dcm(λ + π / 2, π / 2 - ϕ, :ZX)
    else
        throw(ArgumentError("$mount is not a supported topocentric mounting type."))
    end
    pid = axes_id(frames, parent)
    return add_axes_fixedoffset!(frames, name, id, pid, dcm)
end

@fastmath function _geod2pos(h::Number, λ::Number, ϕ::Number, R::Number, f::Number)
    # Get eccentricity from flattening 
    e² = (2 - f) * f

    sϕ, cϕ = sincos(ϕ)
    sλ, cλ = sincos(λ)

    d = R / sqrt(1 - e² * sϕ^2)
    c = (d + h) * cϕ
    s = (1 - e²) * d

    return SA[c*cλ, c*sλ, (s+h)*sϕ]
end

"""
    add_point_surface!(frames, name::Symbol, pointid::Int, parent, axes, 
        λ::Number, ϕ::Number, R::Number, f::Number=0.0, h::Number=0.0)

Add `point` to `frames` as a fixed point on the surface of the `parent` point body. 
The relative  position is specified by the longitude `λ`, the geodetic latitude `ϕ`, 
the reference radius of the ellipsoid `R` and its flattening `f`. 
The altitude over the reference surface of the  ellipsoid `h` defaults to 0. 

!!! warning 
    Axes used here must be a set of body-fixed reference axes for the body represented by 
    `parentid`. This is under user resposibility. 
"""
function add_point_surface!(
    frames::FrameSystem, name::Symbol, id::Int, parent, axes,
    λ::Number, ϕ::Number, R::Number, f::Number=0.0, h::Number=0.0,
)
    pos = _geod2pos(h, λ, ϕ, R, f)
    pid = point_id(frames, parent)
    axid = axes_id(frames, axes)
    return add_point_fixedoffset!(frames, name, id, pid, axid, pos)
end