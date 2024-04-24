function add_axes_topocentric!(
    frames::FrameSystem, name::Symbol, id::Int, parentid::Int, 
    λ::Number, ϕ::Number, type::Symbol
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

    return add_axes_fixedoffset!(frames, name, id, parentid, dcm)
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
    frames::FrameSystem, name::Symbol, id::Int, parentid::Int, axesid::Int, 
    λ::Number, ϕ::Number, R::Number, f::Number=0.0, h::Number=0.0,
)
    pos = _geod2pos(h, λ, ϕ, R, f)
    return add_point_fixed!(frames, name, id, parentid, axesid, pos)
end