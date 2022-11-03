
""" 
    polar_motion(xₚ::N, yₚ::N, t::N)

Polar Motion IAU-2006/2000, CIO Based

### Inputs

- `xp`, `yp` -- Coordinates in radians of the Celestial Intermediate Pole 
                (CIP), with respect to the International Terrestrial Reference
                Frame (ITRF).
                
- `t`  -- Terrestrial Time 

- `sp` -- The Terrestrial Intermediate Origin (TIO) locator, in radians. It provides 
         the position of the TIO on the equatior fo the CIP. 

### Output
Rotation matrix from ITRF to TIRS
"""
function polar_motion(xₚ::Number, yₚ::Number, t::Number)
    sp = tio_locator(t)
    polar_motion(xₚ, yₚ, t, sp)
end

function polar_motion(xₚ::Number, yₚ::Number, ::Number, sp::Number)
    angle_to_dcm(-sp, xₚ, yₚ, :ZYX)
end



"""
    tio_locator(t::N)

Compute the TIO locator `s'`, positioning the Terrestrial Intermediate Origin on 
the equator of the Celestial Intermediate Pole (CIP).

This function approximates the unpredictable motion of the TIO locator s' with 
its secular drift of ~0.47 μas/century. 

### Input
- `tt` -- Terrestrial time 

### Output 
TIO locator at date t
"""
function tio_locator(t::Number)
    return -47e-6*t |> arcsec2rad; # arcseconds
end


""" 
    era_rotm(Tᵤ::N)

Earth Rotation IAU-2006/2000, CIO Based

### Input 
- `Tᵤ` -- Julian UT1 date - 2451545.0

### Output
Rotation matrix from TIRS to CIRS
"""

function era_rotm(Tᵤ::Number) 
    ERA = earth_rotation_angle(Tᵤ)
    angle_to_dcm(-ERA, :Z)
end


"""
    earth_rotation_angle(t::N) where {N <: Number}

Compute the Earth Rotation Angle (ERA), i.e., the angle between the Celestial Intermediate Origin (CIO) 
and the Terrestrial Intermediate Origin (TIO). It uses the fractional UT1 date to gain 
additional precision in the computations (0.002737.. instead of 1.002737..)

### Input 
- `Tᵤ` -- Julian UT1 date - 2451545.0

### Output 
Earth Rotation Angle (ERA) in radians at time Tᵤ
"""
@inline function earth_rotation_angle(Tᵤ::Number)
    return mod2pi(2π * (Tᵤ % 1 + 0.7790572732640 + 0.00273781191135448Tᵤ))
end


"""
    fw2xy(ϵ::Number, ψ::Number, γ::Number, φ::Number)

Compute CIP X and Y coordinates from Fukushima-Williams bias-precession-nutation 
angles.

### Inputs 
- `ϵ` -- F-W angle with IAU 2000A/B nutation corrections. 
- `ψ` -- F-W angle with IAU 2000A/B nutation corrections.
- `γ` -- F-W angle  
- `ϕ` -- F-W angle  

### Outpus 
- `X`, `Y` -- CIP coordinates X, Y

"""
function fw2xy(ϵ::Number, ψ::Number, γ::Number, φ::Number)
    sϵ, cϵ = sincos(ϵ)
    sψ, cψ = sincos(ψ)
    sγ, cγ = sincos(γ)
    sϕ, cϕ = sincos(φ) 

    a = (sϵ*cψ*cϕ - cϵ*sϕ)

    X = sϵ*sψ*cγ - a*sγ
    Y = sϵ*sψ*sγ + a*cγ

    return X, Y
end

"""
    cip_coords
Computes CIP X, Y coordinates 
"""
function cip_coords(m::IAU2006Model, t::Number)
 
    # Computes Fukushima-Williams angles
    γ, ϕ, ψ, ϵ = fw_angles(m, t)

    # Computes IAU 2000 nutation components 
    Δψ, Δϵ = orient_nutation(m, t) 
    
    # Applies IAU-2006 compatible nutations 
    ψₙ = ψ + Δψ
    ϵₙ = ϵ + Δϵ
    
    # Retrieves CIP coordinates 
    fw2xy(ϵₙ, ψₙ, γ, ϕ)
end 


include("constants/cio_locator.jl")

cio_locator(::IAU2006Model, ::Number, ::FundamentalArguments) = ()
build_cio_series(:cio_locator, :IAU2006Model, COEFFS_CIO2006_SP, COEFFS_CIO2006_S)

"""
    cio_locator()

Compute CIO Locator  

Notes: some of the values are slighly different than SOFA but equal to IERS
"""
function cio_locator(m::IAU2006Model, t::Number, x::Number, y::Number)
    fa = FundamentalArguments(t, iau2006a)
    cio_locator(m, t, fa)/1e6*π/648000 - x*y/2
end


function cip_motion(m::IAU2006Model, t::Number)

    # Computes CIP X, Y coordinates 
    x, y = cip_coords(m, t)

    # Computes cio locator `s`
    s = cio_locator(m, t, x, y)
    
    # Retrieves spherical angles E and d. 
    r2 = sqrt(x^2 + y^2)
    E = (r2 > 0) ? atan(y, x) : 0.0 
    d = atan(sqrt(r2 / (1.0 - r2)))

    # This formulation (compared to simplified versions) ensures 
    # that the resulting matrix is orthonormal.
    angle_to_dcm(-E, -d, E+s, :ZYZ)
end


function itrs2gcrs(m::IAU2006Model, t::Number, xₚ::Number, yₚ::Number)

    # Polar motion matrix
    W = polar_motion(xₚ, yₚ, t)

    # Earth Rotation Angle matrix 
    R = era_rotm(t) 

    # Precession-Nutation matrix arising from the CIP motion
    Q = cip_motion(m, t) 

    return Q*R*W 
end