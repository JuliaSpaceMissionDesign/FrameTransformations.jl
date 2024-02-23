export LightTime, PlanetaryAberration

""" 
    light_speed 

Official light speed constant value used in CSPICE 
"""

const light_speed = 299792.458;

# Abstract types definition
abstract type AbstractLightTimeCorrection end
struct LightTimeCorrection <: AbstractLightTimeCorrection end
struct PlanetaryAberrationCorrection <: AbstractLightTimeCorrection end

""" 
    LightTime 

The singleton instance of type `LightTimeCorrection`, used to apply light-time (planetary 
aberration) corrections when computing vectors from the `FrameSystem`.

### See also 
See also [`vector3`](@ref) and [`vector6`](@ref).
"""
const LightTime = LightTimeCorrection()

""" 
    PlanetaryAberration 

The singleton instance of type `PlanetaryAberrationCorrection`, used to apply one-way 
light-time and stellar aberration corrections when computing vectors from the `FrameSystem`.

### See also 
See also [`vector3`](@ref) and [`vector6`](@ref).
"""
const PlanetaryAberration = PlanetaryAberrationCorrection()

# Structure to store and collect useful infos required by 
# light-time corrections routines
struct LTProperties
    aid::Int            # ID of the requested output axes
    ap::Int             # Center point associated to the requested output axes
    dir::Int            # Direction = -1 for Reception, +1 for Transmission
    maxiters::Int
end

# -------------------------------------
# LIGHT TIME 
# -------------------------------------

# Compute position from observer to target corrected for one-way light-time
function light_time_corr3(
    frames::FrameSystem,
    ::LightTimeCorrection,
    ltp::LTProperties,
    from::Int,
    to::Int,
    t::Number,
)

    # Get observer position wrt to SSB 
    pₒ = vector3(frames, 0, from, AXESID_ICRF, t)

    # Compute one-way LT position correction and light-time
    pos, lt = light_time(frames, ltp, to, pₒ, t)

    # Rotate to desired frame
    return _lt_rotation(frames, ltp, from, to, pos, pₒ, t, lt)
end

# Compute state from observer to target corrected for one-way light-time
function light_time_corr6(
    frames::FrameSystem,
    ::LightTimeCorrection,
    ltp::LTProperties,
    from::Int,
    to::Int,
    t::Number,
)

    # Get observer state wrt to SSB 
    xₒ = vector6(frames, 0, from, AXESID_ICRF, t)
    pₒ, vₒ = _posvel(xₒ)

    # Compute one-way LT state correction, light-time and its derivative
    state, lt, dlt = light_time(frames, ltp, to, pₒ, vₒ, t)

    # Rotate to desired frame
    return _lt_rotation(frames, ltp, from, to, state, pₒ, vₒ, t, lt, dlt)
end

# Retrieves the relative position and light-time between a given observer position and a target
function light_time(
    frames::FrameSystem, ltp::LTProperties, to::Int, pobs::AbstractVector, t::Number
)

    # Compute first Light-time guess 
    pₜ = vector3(frames, 0, to, AXESID_ICRF, t)
    lt = light_time(pₜ, pobs)

    iters = 1
    while iters ≤ ltp.maxiters
        # Compute target position at light-time corrected epoch
        pₜ = vector3(frames, 0, to, AXESID_ICRF, t + ltp.dir * lt)

        # Light-time computation 
        lt = light_time(pₜ, pobs)
        iters += 1
    end

    return pₜ - pobs, lt
end

# Compute relative state, light time and its derivative wrt to a given observer state
function light_time(
    frames::FrameSystem,
    ltp::LTProperties,
    to::Int,
    pobs::AbstractVector,
    vobs::AbstractVector,
    t::Number,
)

    # Compute first Light-time guess 
    xₜ = vector6(frames, 0, to, AXESID_ICRF, t)
    pₜ, vₜ = _posvel(xₜ)

    lt = light_time(pₜ, pobs)

    iters = 1
    while iters ≤ ltp.maxiters
        # Compute state position at light-time corrected epoch
        if iters == ltp.maxiters
            xₜ = vector6(frames, 0, to, AXESID_ICRF, t + ltp.dir * lt)
            pₜ, vₜ = _posvel(xₜ)
        else
            # just need the position for intermediate calculations
            pₜ = vector3(frames, 0, to, AXESID_ICRF, t + ltp.dir * lt)
        end

        # Light-time computation 
        lt = light_time(pₜ, pobs)
        iters += 1
    end

    Δp = pₜ - pobs
    Δv = vₜ - vobs

    # Computing light-time derivative
    A = 1 / (light_speed * sqrt(sum(Δp .^ 2)))
    B = dot(Δp, Δv)
    C = dot(Δp, vₜ)

    # Safe-check ensuring the target speed is not above the speed of light 
    # to avoid possible overflow errors! 
    if ltp.dir * C * A > 0.99999999989999999
        throw(
            ErrorException(
                "[Frames] Congrats, you have just won the Nobel Physics Prize! " *
                "The target velocity is above light-speed. You are either a genius or wrong.",
            ),
        )
    end

    dlt = A * B / (1 - ltp.dir * A * C)
    Δv = vₜ * (1 + ltp.dir * dlt) - vobs

    return vcat(Δp, Δv), lt, dlt
end

@inline light_time(a::AbstractVector, b::AbstractVector) = norm(a - b) / light_speed

# Compute Light-time Position Rotation Matrix
function _lt_rotation(
    frames::FrameSystem,
    ltp::LTProperties,
    from::Int,
    to::Int,
    pos::AbstractVector,
    pobs::AbstractVector,
    t::Number,
    lt::Number,
)
    if ltp.aid != AXESID_ICRF
        if is_timefixed(frames, ltp.aid)
            # Rotation does not depend on time! 
            return rotation3(frames, AXESID_ICRF, ltp.aid, t) * pos
        else
            # Get light-time from observer to axes center point
            ltaxes = _lt_axes(frames, ltp, from, to, pobs, t, lt)

            # Compute rotation matrix 
            return rotation3(frames, AXESID_ICRF, ltp.aid, t + ltp.dir * ltaxes) * pos
        end
    end

    # Does not require any rotations 
    return pos
end

# Compute Light-time Position/Velocity Rotation Matrix
function _lt_rotation(
    frames::FrameSystem,
    ltp::LTProperties,
    from::Int,
    to::Int,
    state::AbstractVector,
    pobs::AbstractVector,
    vobs::AbstractVector,
    t::Number,
    lt::Number,
    dlt::Number,
)
    if ltp.aid != AXESID_ICRF
        if is_timefixed(frames, ltp.aid)
            return rotation6(frames, AXESID_ICRF, ltp.aid, t) * state
        else
            # Get light-time from observer to axes center point
            lta, dlta = _lt_axes(frames, ltp, from, to, pobs, vobs, t, lt, dlt)

            # Compute rotation matrix 
            R = rotation6(frames, AXESID_ICRF, ltp.aid, t + ltp.dir * lta)

            # Correct the velocity portion of the rotation matrix to account for 
            # the light-time derivative contribution! 
            @inbounds return Rotation((R[1], (1 + ltp.dir * dlta) * R[2])) * state
        end
    end

    # Does not require any rotations 
    return state
end

# Compute the light-time from observer to axes center point
function _lt_axes(
    frames::FrameSystem,
    ltp::LTProperties,
    from::Int,
    to::Int,
    pobs::AbstractVector,
    t::Number,
    lt::T,
) where {T}
    if ltp.ap == from
        return T(0)
    elseif ltp.ap == to
        return lt
    else
        _, lt = light_time(frames, ltp, ltp.ap, pobs, t)
        return lt
    end
end

# Compute the light-time and its derivative from observer to axes center point
function _lt_axes(
    frames::FrameSystem,
    ltp::LTProperties,
    from::Int,
    to::Int,
    pobs::AbstractVector,
    vobs::AbstractVector,
    t::Number,
    lt::T,
    dlt::T,
) where {T}
    if ltp.ap == from
        return T(0), T(0)
    elseif ltp.ap == to
        return lt, dlt
    else
        _, lt, dlt = light_time(frames, ltp, ltp.ap, pobs, vobs, t)
        return lt, dlt
    end
end

# -------------------------------------
# PLANETARY ABERRATION
# -------------------------------------

# Compute position from observer to target corrected for one-way light-time and 
# stellar aberration
function light_time_corr3(
    frames::FrameSystem,
    ::PlanetaryAberrationCorrection,
    ltp::LTProperties,
    from::Int,
    to::Int,
    t::Number,
)

    # Compute observer state with respect to SSB 
    xₒ = vector6(frames, 0, from, AXESID_ICRF, t)
    pₒ, vₒ = _posvel(xₒ)

    # Correct position for Light-time
    pₜ, lt = light_time(frames, ltp, to, pₒ, t)

    # Compute stellar aberration position contribution
    δp = _stellar_aberration_correction(pₜ, -ltp.dir * vₒ)

    # Rotate to desired frame
    return _lt_rotation(frames, ltp, from, to, pₜ + δp, pₒ, t, lt)
end

# Compute state from observer to target corrected for one-way light-time and 
# stellar aberration
function light_time_corr6(
    frames::FrameSystem,
    ::PlanetaryAberrationCorrection,
    ltp::LTProperties,
    from::Int,
    to::Int,
    t::Number,
)

    # Compute observer state with respect to SSB 
    xₒ = vector9(frames, 0, from, AXESID_ICRF, t)
    pₒ, vₒ, aₒ = _posvelacc(xₒ)

    # Correct state for Light-time
    xₜ, lt, dlt = light_time(frames, ltp, to, pₒ, vₒ, t)

    # Compute stellar aberration state contribution
    δx = _stellar_aberration_correction(xₜ, -ltp.dir * vₒ, -ltp.dir * aₒ)

    # Rotate to desired frame
    return _lt_rotation(frames, ltp, from, to, xₜ + δx, pₒ, vₒ, t, lt, dlt)
end

# Computes the position correction for stellar aberration
function _stellar_aberration_correction(prel::AbstractVector, vobs::AbstractVector)

    # This depends only upon the observer velocity 
    # with respect to the solar system barycenter 

    dᵣ = norm(prel)
    uᵣ = prel / dᵣ

    # Compute sine and cosine of aberration angle! 
    vₚ = vobs - dot(vobs, uᵣ) * uᵣ
    cₐ, sₐ = _aberration_angle(vₚ)

    uᵥ = vₚ / norm(vₚ)
    return ((cₐ - 1) * uᵣ + sₐ * uᵥ) * dᵣ
end

# Computes the state correction for stellar aberration
function _stellar_aberration_correction(
    xrel::AbstractVector, vobs::AbstractVector, aobs::AbstractVector
)

    # This depends only upon the observer velocity and acceleration
    # with respect to the solar system barycenter

    pos, vel = _posvel(xrel)

    dᵣ = norm(pos)
    uᵣ = pos / dᵣ

    δdᵣ = dot(vel, uᵣ)
    δuᵣ = δunitvec(xrel)
    
    # Compute sine and cosine of aberration angle! 
    B = dot(vobs, uᵣ)

    vₚ = vobs - B * uᵣ
    cₐ, sₐ = _aberration_angle(vₚ)

    uᵥ = vₚ / norm(vₚ)

    # Compute aberration position correction
    A = ((cₐ - 1) * uᵣ + sₐ * uᵥ)
    δp = A * dᵣ

    δvₚ = aobs - (dot(aobs, uᵣ) + dot(vobs, δuᵣ)) * uᵣ - B * δuᵣ
    δuᵥ = δunitvec(vcat(vₚ, δvₚ))

    # Compute derivative of aberration angle
    δϕ = 1 / (light_speed * cₐ) * dot(δvₚ, uᵥ)

    # Compute aberration velocity correction
    δv = ((cₐ - 1) * δuᵣ - δϕ * sₐ * uᵣ + sₐ * δuᵥ + δϕ * cₐ * uᵥ) * dᵣ
    δv += A * δdᵣ

    return vcat(δp, δv)
end

# compute sine and cosine of aberration angle, given observer velocity 
function _aberration_angle(vₚ::AbstractVector)
    # Compute sine and cosine of aberration angle! 
    sₐ = norm(vₚ) / light_speed
    cₐ = sqrt(max(0.0, 1 - sₐ^2))

    if cₐ == 0.0
        throw(ErrorException("[Frames] Cosine of aberration angle cannot be <= 0.0"))
    end

    return cₐ, sₐ
end

# Utilities functions to split state vector
@inline _posvel(x::AbstractVector) = @inbounds SA[x[1], x[2], x[3]], SA[x[4], x[5], x[6]]

@inline function _posvelacc(x::AbstractVector)
    @inbounds SA[x[1], x[2], x[3]], SA[x[4], x[5], x[6]], SA[x[7], x[8], x[9]]
end
