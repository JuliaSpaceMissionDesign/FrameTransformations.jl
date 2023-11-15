export orient_gmst, orient_gast

# Build the function to compute the EE complementary terms
include("constants/eq_origins.jl");
build_series(:ee_complementary, :IAU2000Model, [COEFFS_EECT])


"""
    orient_gmst(::IAUModel, t::Number)

Compute the Greenwich Mean Sidereal Time (GMST), in radians, according to the IAU Model `m`, 
given time `tt` expressed in `TT` Julian centuries since `J2000`. 

!!! note 
    This function is valid for the IAU2000 and IAU2006 models only.

!!! note 
    This function computes the Earth Rotation Angle (ERA) by automatically converting TT to 
    UT1. Lower-level interfaces are also available to minimise the number of repeated 
    computations. 

### See also 
See also [`orient_gast`](@ref).
"""
function orient_gmst(m::IAUModel, tt::Number)

    # Transform TT centuries to UT1 days
    ut1_d = iers_tt_to_ut1(tt * Tempo.CENTURY2DAY)

    # Compute the Earth Rotation Angle 
    θ = earth_rotation_angle(ut1_d)

    # Compute the GMST
    return orient_gmst(m, tt, θ)

end

"""
    orient_gmst(m::IAU1980Model, ut1::Number)

Compute the Greenwich Mean Sidereal Time (GMST), in radians, according to the IAU 1980 
models, given time `ut1` expressed in `UT1` Julian centuries since `J2000`. 

### References 
- [ERFA gmst82](https://github.com/liberfa/erfa/blob/master/src/gmst82.c) routine.
"""
function orient_gmst(::IAU1980Model, ut1::Number)

    # Coefficients for the IAU 1982 GMST-UT1 model. The first component has been adjusted 
    # of 12 hours because UT1 starts at noon.
    A = -19089.45159
    B = 8640184.812866
    C = 0.093104
    D = -6.2e-6

    # Fractional part of UT1, in seconds
    f = Tempo.DAY2SEC*mod(ut1*Tempo.CENTURY2DAY, 1)

    # seconds to radians conversion 
    s2r = 7.272205216643039903848712e-5

    # Compute GMST
    return mod2pi(s2r*(@evalpoly(ut1, A, B, C, D) + f))

end

"""
    orient_gmst(m::IAUModel, tt::Number, θ::Number)

Compute the Greenwich Mean Sidereal Time (GMST), in radians, according to the IAU Model `m`, 
given the Earth Rotation Angle (ERA) `θ`, in radians and the time `tt` expressed in `TT` 
Julian centuries since `J2000`. 

The function has been implemented for the IAU2000 and IAU2006 models.

### References 
- [ERFA gmst00](https://github.com/liberfa/erfa/blob/master/src/gmst00.c) routine.
- [ERFA gmst06](https://github.com/liberfa/erfa/blob/master/src/gmst06.c) routine.
"""
function orient_gmst(::IAU2006Model, tt::Number, θ::Number)

    # Evaluate interpolating series
    p = @evalpoly(
        tt, 
        0.014506, 
     4612.156534,   
        1.3915817, 
       -0.00000044,
       -0.000029956, 
       -0.0000000368
    )

    # Compute GMST 
    return mod2pi(θ + arcsec2rad(p))

end

function orient_gmst(::IAU2000Model, tt::Number, θ::Number)

    # Evaluate interpolating series
    p = @evalpoly(
        tt, 
        0.014506, 
     4612.15739966,   
        1.39667721, 
       -0.00009344,
        0.00001882
    )

    # Compute GMST 
    return mod2pi(θ + arcsec2rad(p))

end


"""
    orient_gast(m::IAUModel, t::Number)

Compute the Greenwich Apparent Sidereal Time (GAST), in radians, given time `t` as `TT`
Julian  centuries since `J2000` according to the IAU Model `m`.

!!! note 
    For the IAU2000B model, as an approximation ERFA uses UT1 instead of TDB (or TT) to 
    compute the precession component of GMST and the equation of the equinoxes. This 
    approximation is not performed in this framework.

!!! note 
    This function computes the Earth Rotation Angle (ERA) by automatically converting TT to 
    UT1. Lower-level interfaces are also available to minimise the number of repeated 
    computations. 

### See also 
See also [`orient_gmst`](@ref).
"""
function orient_gast(m::IAUModel, t::Number)

    # Transform TT centuries to UT1 days
    ut1_d = iers_tt_to_ut1(t * Tempo.CENTURY2DAY)

    # Compute the Earth Rotation Angle 
    θ = earth_rotation_angle(ut1_d)

    # Compute GAST 
    return orient_gast(m, t, θ)
    
end

"""
    orient_gast(m::IAUModel, t::Number, θ::Number)

Compute the Greenwich Apparent Sidereal Time (GAST), in radians, given time `t` as `TT` 
Julian centuries since `J2000` and the Earth Rotation Angle (ERA) `θ`, in radians, according 
to the IAU Model `m`.

The function has been implemented for the IAU2000 and IAU2006 models.

### References 
- [ERFA gst00a](https://github.com/liberfa/erfa/blob/master/src/gst00a.c) routine.
- [ERFA gst06a](https://github.com/liberfa/erfa/blob/master/src/gst06a.c) routine.
"""
function orient_gast(m::IAU2006Model, t::Number, θ::Number)

    # Compute the equations of the origins
    eors = origins_equation(m, t)  
    return mod2pi(θ - eors)

end

function orient_gast(m::IAU2000Model, t::Number, θ::Number)

    # Compute the Greenwich Mean Sidereal time 
    gmst = orient_gmst(m, t, θ)

    # Compute the equations of the equinoxes
    ee = equinoxes_equation(m, t)

    # Comput GAST, in radians
    return mod2pi(gmst + ee)

end

# function orient_gast(::IAU1980Model)

#     # Compute the Greenwich Mean Sidereal time 
#     gmst = orient_gmst(m, t)

#     # Compute the equations of the equinoxes
#     ee = equinoxes_equation(m, t)

#     # Comput GAST, in radians
#     return mod2pi(GMST + ee)

# end


"""
    equinoxes_equation(m::IAUModel, tt::Number)

Compute the Equation of the Equinoxes, in radians, according to the IAU Model `m` given 
time `tt` expressed in `TT` Julian centuries since J2000. 

This function has been implemented for the IAU2006 and IAU2000 models.

!!! note 
    This function neglects the difference between `TT` and `TDB`.

### References 
- ERFA [ee00a](https://github.com/liberfa/erfa/blob/master/src/ee00a.c) routine.
- ERFA [ee06a](https://github.com/liberfa/erfa/blob/master/src/ee06a.c) routine.
"""
function equinoxes_equation(m::IAU2006Model, tt::Number)

    # Compute the Earth Rotation Angle (ERA) 
    θ = earth_rotation_angle(-Tempo.DJ2000)

    # Compute GMST and GAST, in radians 
    gast = orient_gast(m, tt, θ)
    gmst = orient_gmst(m, tt, θ)

    # Equation of the equinoxes
    return mod2pi(gast - gmst)

end

# Time expressed in TT julian centuries since J2000
function equinoxes_equation(m::IAU2000Model, tt::Number)

    # we neglect the difference between TT and TDB for the FAs...

    # Compute precession-rate adjustments 
    _, Δϵₚ = precession_rate(m, tt)

    # Compute the mean obliquity 
    ϵₐ = orient_obliquity(iau1980, tt) + Δϵₚ

    # Nutation in logitutude 
    Δψ, _ = orient_nutation(m, tt)

    # Compute the Fundamental Arguments using the IAU2000A model because we need 
    # the associated expressions for the Luni-solar arguments
    fa = FundamentalArguments(tt, iau2000a)

    # Equation of the equinoxes
    return Δψ*cos(ϵₐ) + ee_complementary(m, tt, fa)

end

# # Time expressed in TDB julian centuries since J2000
# function equinoxes_equation(m::IAU1980Model, t::Number)

#     # Longitude of the mean ascending node of the lunar orbit on the ecliptic, 
#     # measured from the mean equinox of date, in radians
#     ω = mod2pi(arcsec2rad(@evalpoly(t, -482890.539, 7.455, 0.008)) + 2π*mod(-5t, 1))

#     # Compute the nutation components in longitude and obliquity
#     # TODO: missing orient_nutation per il 1980 missa 

#     # Compute the mean obliquity 
#     ϵₐ = orient_obliquity(m, t)

#     # Equations of the equinoxes 
#     return Δψ*cos(ϵₐ) + arcsec2rad(0.00264*sin(ω) + 0.000063*sin(2ω))

# end

"""
    origins_equation(m::IAU2006Model, t::Number)

Compute the Equation of the Origins (EO), in radians, following the IAU2006 precession and 
IAU2000A nutation models, given time `t` expressed in `TT` Julian centuries since J2000.0. 

### References
- Wallace P. T. and Capitaine N. (2006), Precession-nutation procedures consistent with 
  IAU 2006 resolutions, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
- ERFA [eo06a](https://github.com/liberfa/erfa/blob/master/src/eo06a.c) routine.
"""
function origins_equation(m::IAU2006Model, t::Number)

    # Compute the bias-precession-nutation matrix 
    bpn = orient_bias_precession_nutation(m, t)

    # Compute the CIP coordinates 
    x, y = bpn2xy(bpn)    

    # Get the CIO locator 
    s = cio_locator(m, t, x, y)

    @inbounds begin 
        ax = x / (1 + bpn[3, 3])
        xs = 1 - ax*x 
        ys = - ax * bpn[3, 2]
        zs = -x 

        p = bpn[1, 1]*xs + bpn[1, 2]*ys + bpn[1, 3] * zs 
        q = bpn[2, 1]*xs + bpn[2, 2]*ys + bpn[2, 3] * zs
        
    end

    return ((p != 0) || (q != 0)) ? s - atan(q, p) : s
        
end

