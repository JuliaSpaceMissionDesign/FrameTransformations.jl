export orient_rot3_itrf_to_gcrf

# Nominal Earth angular velocity  
const ωₑ = 7.292_115_146_706_979e-5


"""
    earth_rotation_rate()
    
Compute the nominal Earth angular velocity. 
"""
earth_rotation_rate() = earth_rotation_rate(0.0)


""" 
    earth_rotation_rate(LOD::Number)

Compute the true angular velocity of the Earth accounting for the Length of the Day, i.e., 
the instantaneous rate of change of UT1 with respect to a uniform time scale. 
"""
function earth_rotation_rate(LOD::Number)
    return ωₑ*(1-LOD/86400)
end


"""
    tio_locator(t::Number)

Compute the TIO locator `s'` at date, positioning the Terrestrial Intermediate Origin on 
the equator of the Celestial Intermediate Pole (CIP) at time `t` expressed as `TT` Julian 
centuries since J2000. 

This function approximates the unpredictable motion of the TIO locator s' with its secular 
drift of ~0.47 μas/century. 

### References
- Luzum, B. and Petit G. (2012), The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- Lambert, S. and Bizouard C. (2002), Positioning the Terrestrial Ephemeris Origin in the 
  Terrestrial Reference Frame, [DOI: 10.1051/0004-6361:20021139](https://www.aanda.org/articles/aa/pdf/2002/40/aa2747.pdf)
"""
function tio_locator(t::Number)
    return -47e-6*t |> arcsec2rad; # arcseconds
end


""" 
    polar_motion(t::Number, xₚ::Number, yₚ::Number)

Compute the Polar Motion rotation matrix from ITRF to TIRS, according to the 
IERS 2010 Conventions, at time `t` expressed in `TT` Julian centuries since [`J2000`](@ref). 
The function requires `xp` and `yp`, the Celestial Intermediate Pole (CIP) coordinates with 
respect to the International Celestial Reference Frame (ITFR).

### References 
- Luzum, B. and Petit G. (2012), The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
"""
function polar_motion(t::Number, xₚ::Number, yₚ::Number,)
    sp = tio_locator(t)
    angle_to_dcm(yₚ, xₚ, -sp, :XYZ)
end


"""
    earth_rotation_angle(t::Number)

Compute the Earth Rotation Angle (ERA) in radians, i.e., the angle between the Celestial 
Intermediate Origin (CIO) and the Terrestrial Intermediate Origin (TIO) at time `t` 
expressed as UT1 days since [`J2000`](@ref).

!!! note 
    The function uses the fractional UT1 date to gain additional precision in the 
    computations (0.002737.. instead of 1.002737..)

### References 
- Luzum, B. and Petit G. (2012), The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/era00.c) library
"""
@inline function earth_rotation_angle(t::Number)
    f = t % 1.0
    return mod2pi(2π * (f + 0.7790572732640 + 0.00273781191135448t))
end

""" 
    era_rotm(t::Number)

Compute the TIRS to CIRS Earth Rotation matrix, according to the IERS 2010 
conventions at time `t` expressed as UT1 days since [`J2000`](@ref).

### References
- Luzum, B. and Petit G. (2012), The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
"""
function era_rotm(t::Number) 
    angle_to_dcm(-earth_rotation_angle(t), :Z)
end


"""
    fw2xy(ϵ::Number, ψ::Number, γ::Number, φ::Number)

Compute the CIP X and Y coordinates from Fukushima-Williams bias-precession-nutation 
angles, in radians.

### Inputs 
- `ϵ` -- F-W angle with IAU 2006A/B nutation corrections. 
- `ψ` -- F-W angle with IAU 2006A/B nutation corrections.
- `γ` -- F-W angle  
- `ϕ` -- F-W angle

### References
- Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent with 
  IAU 2006 resolutions, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
"""
function fw2xy(γ::Number, φ::Number, ψ::Number, ϵ::Number)
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
    bpn2xy(A::AbstractMatrix)

Compute the CIP X and Y coordinates from the bias-precession-nutation matrix, in radians.

### References 
- Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent with 
  IAU 2006 resolutions, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
"""
bpn2xy(A::AbstractMatrix) = @inbounds A[3, 1], A[3, 2]

# Generate cip_coords function for the simplified CPNC model! 
include("constants/cip_cpnc.jl")
build_cio_series(
    :cip_coords, :CPNC, 
    [COEFFS_CPNC_XP, COEFFS_CPNC_YP], [COEFFS_CPNC_X, COEFFS_CPNC_Y]
)


"""
    cip_coords(m::IAUModel, t::Number)

Computes the CIP X, Y coordinates, in radians, according to the IAU model `m` at time `t` 
expressed in `TT` Julian Centuries since [`J2000`](@ref).

This function has been implemented for the `IAU2000`, `IAU2006` and the `CPN` models.

### References
- Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent with 
  IAU 2006 resolutions, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
- Capitaine N. and Wallace P. T. (2008), Concise CIO based precession-nutation formulations
"""
function cip_coords(m::IAU2006Model, t::Number)
 
    # Computes Fukushima-Williams angles
    γ, ϕ, ψ, ϵ = fw_angles(m, t)

    # Computes IAU 2000 nutation components 
    Δψ, Δϵ = orient_nutation(m, t) 

    # Retrieves CIP coordinates by applying IAU-2006 compatible nutations 
    fw2xy(γ, ϕ, ψ + Δψ, ϵ + Δϵ)
end 

function cip_coords(m::IAU2000Model, t::Number)
    # Extract CIP coordinates from the IAU-2000 bias-precession-nutation matrix
    bpn2xy(orient_bias_precession_nutation(m, t))
end 

function cip_coords(::CPND, t::Number)

    μas2rad = 1e-6*π/648000

    # Approximated fundamental arguments as linear function of time 
    Ω = 2.182439196616 - 33.7570459536t
    A = -2.776244621014 + 1256.6639307381t

    sΩ, cΩ = sincos(Ω)
    sA, cA = sincos(A)

    X = @evalpoly(t, 0.0, 2004191898.0, -429782.9, -198618.34)
    X -= 6844318.0sΩ + 523908.0sA

    Y = @evalpoly(t, 0.0, 0.0, -22407275.0)
    Y += 9205236.0cΩ + 573033.0cA

    X *= μas2rad
    Y *= μas2rad

    return X, Y
end

function cip_coords(m::CPNC, t::N) where {N <: Number}
    
    # Computes only Luni-Solar Fundamental Arguments 

    z = N(0) 

    # Computes only Luni-Solar Fundamental Arguments 
    fa = FundamentalArguments(
          LuniSolarArguments(t, m)..., z, z, z, z, z, z, z, z, z
    )


    cip_coords(m, t, fa)
end



include("constants/cio_locator00.jl")
include("constants/cio_locator06.jl")

# cio_locator(::IAUModel, ::Number, ::FundamentalArguments) = ()
build_cio_series(:cio_locator, :IAU2000Model, [COEFFS_CIO2000_SP], [COEFFS_CIO2000_S])
build_cio_series(:cio_locator, :IAU2006Model, [COEFFS_CIO2006_SP], [COEFFS_CIO2006_S])

"""
    cio_locator(m::IAUModel, t::Number, x::Number, y::Number)

Compute the CIO Locator `s` in radians, according to the IAU Model `m`, given the CIP 
coordinates `X` and `Y` at time `t` expressed in `TT` Julian centuries since [`J2000`](@ref)

The function has been implemented for the `IAU2000`, `IAU2006` and the `CPN` models.

### References 
- Luzum, B. and Petit G. (2012), The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent with 
  IAU 2006 resolutions, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
- Capitaine N. and Wallace P. T. (2008), Concise CIO based precession-nutation formulations
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/s00.c) library
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/s06.c) library
"""
function cio_locator(m::IAUModel, t::Number, x::Number, y::Number)
    fa = FundamentalArguments(t, iau2006a)
    cio_locator(m, t, fa) - x*y/2
end

function cio_locator(::CPNC, t::Number, x::Number, y::Number)
    # Simplified model! 
    Ω = 2.182439196616 - 33.7570459536t

    s = @evalpoly(t, 0.0, 3809, 0.0, -72574.0)
    s -= 2641*sin(Ω)

    # Transform s from μas to radians 
    s *=  1e-6*π/648000

    return s - x*y/2

end

@inline cio_locator(::CPND, t::Number, x::Number, y::Number) = 0.0


""" 
    xys2m(x::Number, y::Number, s::Number)

Compute the Intermediate-to-Celestial matrix given the CIP `x`, `y' coordinates and the CIO 
locator `s`, all in radians.

### References
- Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent with 
  IAU 2006 resolutions, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/c2ixys.c) library

"""
function xys2m(x::Number, y::Number, s::Number)
    # Retrieves spherical angles E and d. 
    r2 = x^2 + y^2
    E = (r2 > 0) ? atan(y, x) : 0.0 
    d = atan(sqrt(r2 / (1.0 - r2)))

    # This formulation (compared to simplified versions) ensures 
    # that the resulting matrix is orthonormal.
    angle_to_dcm(E+s, -d, -E, :ZYZ)
end


"""
    cip_motion(m::IAUModel, t::Number, dx::Number=0.0, dy::Number=0.0)

Compute the CIRS to GCRS rotation matrix, according to the IAU Model `m`, at time `t`
expressed in `TT` Julian centuries since [`J2000`](@ref). Optional IERS corrections for 
free-core nutation and time depedented effects can be provided through `dx` and `dy`. 

### References 
- Luzum, B. and Petit G. (2012), The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent with 
  IAU 2006 resolutions, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
"""
function cip_motion(m::IAUModel, t::Number, dx::Number=0.0, dy::Number=0.0)
    # Compute CIP coordinates
    x, y = cip_coords(m, t)

    # Apply free-core nutation corrections
    x += dx 
    y += dy 

    # Compute CIO Locator
    s = cio_locator(m, t, x, y)

    # Form intermediate-to-celestial matrix 
    xys2m(x, y, s)
end


# General functions to dispatch generic order derivatives!
function _dna_itrf_to_gcrf(m::IAUModel, fn::Function, t::Number)

    utc_s = Tempo.apply_offsets(Tempo.TIMESCALES, t, TT, UTC)
    ut1 = Tempo.apply_offsets(Tempo.TIMESCALES, utc_s, UTC, UT1)/Tempo.DAY2SEC
    
    utc_d = utc_s/Tempo.DAY2SEC

    # Compute pole coordinates 
    xₚ = interpolate(IERS_EOP.x, utc_d) |> arcsec2rad
    yₚ = interpolate(IERS_EOP.y, utc_d) |> arcsec2rad

    # Compute dX, dY 
    dX = 1e-3*interpolate(IERS_EOP.dX, utc_d) |> arcsec2rad
    dY = 1e-3*interpolate(IERS_EOP.dY, utc_d) |> arcsec2rad

    # Compute LOD 
    LOD = 1e-3*interpolate(IERS_EOP.LOD, utc_d)

    return fn(m, t, ut1, xₚ, yₚ, dX, dY, LOD)

end

function _dnb_itrf_to_gcrf(m::IAUModel, fn::Function, t::Number)

    # Convert TT secs since J2000 to TT days
    ttd = t/Tempo.DAY2SEC

    # Compute pole coordinates 
    xₚ = interpolate(IERS_EOP.x_TT, ttd) |> arcsec2rad
    yₚ = interpolate(IERS_EOP.y_TT, ttd) |> arcsec2rad

    # Transform UT1 to TT
    offset = interpolate(IERS_EOP.UT1_TT, ttd)
    ut1 = ttd + offset/Tempo.DAY2SEC

    return fn(m, t, ut1, xₚ, yₚ, 0.0, 0.0)
end

function _dnd_itrf_to_gcrf(m::IAUModel, fn::Function, t::Number)

    # Convert TT secs since J2000 to TT days
    ttd = t/Tempo.DAY2SEC

    # Transform UT1 to TT
    offset = interpolate(IERS_EOP.UT1_TT, ttd)
    ut1 = ttd + offset/Tempo.DAY2SEC

    return fn(m, t, ut1, 0.0, 0.0, 0.0, 0.0)
end


"""
    orient_rot3_itrf_to_gcrf(m::IAUModel, t::Number)

Compute the rotation matrix from `ITRF` to `GCRF` at time `t` expressed as 
TT seconds since [`J2000`](@ref), according to the IAU Model `m`, as follows:

- **IAU2000A**: the pole coordinates (xₚ, yₚ) and the free-core nutation and time corrections 
    to the CIP coordinates (dX, dY) are interpolated from the latest released IERS EOP data. 
    The precession-nutation matrix is computed using the full IAU 2000A model.

- **IAU2000B**: only the pole coordinates (xₚ, yₚ) are interpolated from the latest EOP data. 
    The Free Core Nutation (FCN) corrections dX, dY are neglected. The precession-nutation 
    matrix is computed following the IAU 2000 model but with truncated expressions for the 
    nutation corrections. 

- **IAU2006A**: the pole coordinates (xₚ, yₚ) and the free-core nutation and time corrections 
    to the CIP coordinates (dX, dY) are interpolated from the latest released IERS EOP data. 
    The precession-nutation matrix is computed using the full IAU 2006/2000A model.

- **IAU2006B**: only the pole coordinates (xₚ, yₚ) are interpolated from the latest EOP data. 
    The Free Core Nutation (FCN) corrections dX, dY are neglected. The precession-nutation
    matrix is computed following the IAU 2006A model but with truncated expressions for the 
    nutation corrections. 

- **CPNc**: a concise model with a cut-off at 2.5 mas of the X and Y series, delivering a 
    worst-case accuracy of about 15 mas between 1995-2050. It does not take into account the 
    Free Core Nutation (~0.2 mas). 
                    
- **CPNd**: an extremely concise formulation with an accuracy of about 1 arcsec between 1995 
    and 2050. It neglects polar-motion (~0.25 arcsec), the FCN corrections and the CIO locator. 

### References 
- Luzum, B. and Petit G. (2012), The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent with 
  IAU 2006 resolutions, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
- Capitaine N. and Wallace P. T. (2008), Concise CIO based precession-nutation formulations
"""
function orient_rot3_itrf_to_gcrf(m::Union{<:IAU2000A, <:IAU2006A}, t::Number)

    # Find UT1 and UTC dates
    utc_s = Tempo.apply_offsets(Tempo.TIMESCALES, t, TT, UTC)
    ut1 = Tempo.apply_offsets(Tempo.TIMESCALES, utc_s, UTC, UT1)/Tempo.DAY2SEC

    utc_d = utc_s/Tempo.DAY2SEC

    # Compute pole coordinates 
    xₚ = interpolate(IERS_EOP.x, utc_d) |> arcsec2rad
    yₚ = interpolate(IERS_EOP.y, utc_d) |> arcsec2rad

    # Compute dX, dY 
    dX = 1e-3*interpolate(IERS_EOP.dX, utc_d) |> arcsec2rad
    dY = 1e-3*interpolate(IERS_EOP.dY, utc_d) |> arcsec2rad

    return orient_rot3_itrf_to_gcrf(m, t, ut1, xₚ, yₚ, dX, dY)
end


"""
    orient_rot6_itrf_to_gcrf(m::IAUModel, t::Number)

Compute the rotation matrix from `ITRF` to `GCRF` and its derivative at time `t` 
expressed as TT seconds since [`J2000`](@ref), according to the the IAU Model `m`.
"""
orient_rot6_itrf_to_gcrf

"""
    orient_rot9_itrf_to_gcrf(m::IAUModel, t::Number)

Compute the rotation matrix from `ITRF` to `GCRF` and its time derivatives up to order 2 at 
time `t` expressed as TT seconds since [`J2000`](@ref), according to the IAU Model `m`.
"""
orient_rot9_itrf_to_gcrf


"""
    orient_rot12_itrf_to_gcrf(m::IAUModel, t::Number)

Compute the rotation matrix from `ITRF` to `GCRF` and its time derivatives up to order 3 at 
time `t` expressed as TT seconds since [`J2000`](@ref), according to the IAU Model `m`.
"""
orient_rot12_itrf_to_gcrf


for (i, fun) in enumerate([
                :orient_rot3_itrf_to_gcrf, :orient_rot6_itrf_to_gcrf,
                :orient_rot9_itrf_to_gcrf, :orient_rot12_itrf_to_gcrf
                ])

    if i > 1
        @eval begin 
            @inline function ($fun)(m::Union{<:IAU2000A, <:IAU2006A}, t::Number)
                _dna_itrf_to_gcrf(m, $fun, t)
            end
        end
    
    end

    @eval begin 
        @inline function ($fun)(m::Union{<:IAU2000B, <:IAU2006B, <:CPNC}, t::Number)
            _dnb_itrf_to_gcrf(m, $fun, t)
        end
    end

    @eval begin 
        @inline function ($fun)(m::CPND, t::Number)
            _dnd_itrf_to_gcrf(m, $fun, t)
        end
    end

end


"""
    orient_rot3_itrf_to_gcrf(m::IAUModel, tt, ut1, xₚ, yₚ, dX=0.0, dY=0.0)

Compute the rotation matrix from ITRF to GCRF according to the IAU Model `m`, at time `tt` 
and `ut1` expressed in TT seconds and `UT1` days since [`J2000`](@ref), respectively.  

This function has been implemented for `IAU2000`, `IAU2006` and `CPN` models.

!!! note 
    All the input quantities `xₚ`, `yₚ`, `dX` and `dY` must be expressed in radians
"""
function orient_rot3_itrf_to_gcrf(m::IAUModel, tt::Number, ut1::Number, xₚ::Number, 
            yₚ::Number, dX::Number=0.0, dY::Number=0.0)

    # Convert TT since J2000 from seconds to centuries
    tt_c = tt/Tempo.CENTURY2SEC

    W = polar_motion(tt_c, xₚ, yₚ)
    R = era_rotm(ut1)
    Q = cip_motion(m, tt_c, dX, dY)

    return Q*R*W 

end


"""
    orient_rot6_itrf_to_gcrf(m::IAUModel, tt, ut1, xₚ, yₚ, dX=0.0, dY=0.0, LOD=0.0)

Compute the rotation matrix from ITRF to GCRF and its derivative, according to the IAU Model 
`m`, at time `tt` and `ut1` expressed in TT seconds and `UT1` days since [`J2000`](@ref), 
respectively. 

This function has been implemented for `IAU2000`, `IAU2006` and `CPN` models.

!!! note 
    All the input quantities `xₚ`, `yₚ`, `dX` and `dY` must be expressed in radians
"""
function orient_rot6_itrf_to_gcrf(m::IAUModel, tt::Number, ut1::Number , xₚ::Number, 
            yₚ::Number, dX::Number=0.0, dY::Number=0.0, LOD::Number=0.0)

    # Convert TT since J2000 from seconds to centuries
    tt_c = tt/Tempo.CENTURY2SEC

    W = polar_motion(tt_c, xₚ, yₚ)
    R = era_rotm(ut1)
    Q = cip_motion(m, tt_c, dX, dY)

    ωe = SVector(0.0, 0.0, earth_rotation_rate(LOD))
    Ω = skew(ωe)

    QR = Q*R

    D = QR*W
    δD = QR*Ω*W

    return D, δD
end


"""
    orient_rot9_itrf_to_gcrf(m::IAUModel, tt, ut1, xₚ, yₚ, dX=0.0, dY=0.0, LOD=0.0)

Compute the rotation matrix from ITRF to GCRF and its derivatives up to order 2, according 
to the IAU Model `m`, at time `tt` and `ut1` expressed in TT seconds and `UT1` days since 
[`J2000`](@ref), respectively. 

This function has been implemented for `IAU2000`, `IAU2006` and `CPN` models.

!!! note 
    All the input quantities `xₚ`, `yₚ`, `dX` and `dY` must be expressed in radians
"""
function orient_rot9_itrf_to_gcrf(m::IAUModel, tt::Number, ut1::Number, xₚ::Number, 
            yₚ::Number, dX::Number=0.0, dY::Number=0.0, LOD::Number=0.0)

    # Convert TT since J2000 from seconds to centuries
    tt_c = tt/Tempo.CENTURY2SEC

    W = polar_motion(tt_c, xₚ, yₚ)
    R = era_rotm(ut1)
    Q = cip_motion(m, tt_c, dX, dY)

    ωe = SVector(0.0, 0.0, earth_rotation_rate(LOD))
    Ω = skew(ωe)

    QR = Q*R 
    QRΩ = QR*Ω

    D = QR*W
    δD = QRΩ*W
    δ²D = QRΩ*Ω*W

    return D, δD, δ²D
end


"""
    orient_rot12_itrf_to_gcrf(m::IAUModel, tt, ut1, xₚ, yₚ, dX=0.0, dY=0.0, LOD=0.0)

Compute the rotation matrix from ITRF to GCRF and its derivatives up to order 3, according 
to the IAU Model `m`, at time `tt` and `ut1` expressed in TT seconds and `UT1` days since 
[`J2000`](@ref), respectively. 

This function has been implemented for `IAU2000`, `IAU2006` and `CPN` models.

!!! note 
    All the input quantities `xₚ`, `yₚ`, `dX` and `dY` must be expressed in radians
"""
function orient_rot12_itrf_to_gcrf(m::IAUModel, tt::Number, ut1::Number, xₚ::Number, 
            yₚ::Number, dX::Number=0.0, dY::Number=0.0, LOD::Number=0.0)

    # Convert TT since J2000 from seconds to centuries
    tt_c = tt/Tempo.CENTURY2SEC

    W = polar_motion(tt_c, xₚ, yₚ)
    R = era_rotm(ut1)
    Q = cip_motion(tt_c, m, dX, dY)

    ωe = SVector(0.0, 0.0, earth_rotation_rate(LOD))
    Ω = skew(ωe)

    QR = Q*R 
    QRΩ = QR*Ω
    QRΩ² = QRΩ*Ω

    D = QR*W
    δD = QRΩ*W
    δ²D = QRΩ²*W
    δ³D = QRΩ²*Ω*W

    return D, δD, δ²D, δ³D
end