
""" 
    polar_motion(xₚ::N, yₚ::N, s::N)

Polar Motion IAU-2006/2000, CIO Based

### Inputs

- `xp`, `yp` -- Coordinates in radians of the Celestial Intermediate Pole 
                (CIP), with respect to the International Terrestrial Reference
                Frame (ITRF).
                
- `s` -- The Terrestrial Intermediate Origin (TIO) locator, in radians. It provides 
         the position of the TIO on the equatior fo the CIP. 

### Output
Rotation matrix from ITRF to TIRS
"""
function polar_motion(xₚ::Number, yₚ::Number, s::Number)
    sx, cx = sincos(xₚ)
    sy, cy = sincos(yₚ)
    ss, cs = sincos(s)

    # required because an inexact rounding error might be throw if N is an Int 
    SMatrix{3, 3, typeof(sx), 9}(cx*cs, cx*ss, sx, 
                                (-cy*ss + sy*sx*cs), (cy*cs + sy*sx*ss),-sy*cx,  
                                (-sy*ss-cy*sx*cs), (sy*cs - cy*sx*ss), cy*cx)
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
    s, c = sincos(-earth_rotation_angle(Tᵤ) )

    # required because an inexact rounding error might be throw if N is an Int 
    SMatrix{3, 3, typeof(s), 9}(c, -s, 0, s, c, 0, 0, 0, 1)
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
    Y = sϵ*sψ*sγ - a*cγ

    return X, Y
end

"""
    cip_coords
Computes CIP X, Y coordinates 
"""
function cip_coords(m::IAU2006Model, t::Number)
 
    # Computes Fukushima-Williams angles
    ϵ, ψ, γ, ϕ = fw_angles(m, t)

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
    cio_locator(m, t, fa)*1e-6*π/648000 - x*y/2
end


function cirs2gcrs(m::IAU2006Model, t::Number)

    # Computes CIP X, Y coordinates 
    x, y = cip_coords(m, t)

    # Computes cio locator `s`
    s = cio_locator(m, t, x, y) 

    # Computes rotation matrix from CIRS to GCRS
    a = 0.5 + 1/8*(x^2 + y^2)

    Rs = SMatrix{3, 3, typeof(a), 9}(1-a*x^2, -a*x*y,    -x, 
                                     -a*x*y,  1-a*y^2,   -y, 
                                     x, y,   -a*(x^2 + y^2))

    sₛ, cₛ = sincos(s) 

    Rs*SMatrix{3, 3, typeof(a), 9}(cₛ, -sₛ, 0, sₛ, cₛ, 0, 0, 0, 1)

end


function itrs2gcrs(m::IAU2006Model, t::Number, xₚ::Number, yₚ::Number)

    # Polar motion matrix
    sp = tio_locator(t)
    W = polar_motion(xₚ, yₚ, sp)

    # Earth Rotation Angle matrix 
    R = era_rotm(t) 

    # Precession-Nutation matrix 
    Q = cirs2gcrs(m, t) 

    return Q*R*W 
end