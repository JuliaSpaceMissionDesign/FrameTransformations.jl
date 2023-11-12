export gmst


"""
    gmst(m::IAUModel, t::Number, θ::Number)

Compute the Greenwich Mean Sidereal Time (GMST), in radians, according to the IAU Model `m`, 
given the Earth Rotation Angle (ERA) `θ`, in radians and the time `t` expressed in `TT` 
Julian centuries since `J2000`. 

The function has been implemented for the IAU2000 and IAU2006 models.

### References 
- [ERFA gmst00](https://github.com/liberfa/erfa/blob/master/src/gmst00.c) routine.
- [ERFA gmst06](https://github.com/liberfa/erfa/blob/master/src/gmst06.c) routine.
"""
function gmst(::IAU2006Model, t::Number, θ::Number, )

    # Evaluate interpolating series
    p = @evalpoly(
        t, 
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

function gmst(::IAU2000Model, t::Number, θ::Number, )

    # Evaluate interpolating series
    p = @evalpoly(
        t, 
        0.014506, 
     4612.15739966,   
        1.39667721, 
       -0.00009344,
       -0.00001882
    )

    # Compute GMST 
    return mod2pi(θ + arcsec2rad(p))

end

"""
Compute the Greenwich Mean Sidereal Time (GMST), in radians, according to the IAU 1980 
models, given time `t` expressed in `TT` Julian centuries since `J2000`. 

### References 
- [ERFA gmst82](https://github.com/liberfa/erfa/blob/master/src/gmst82.c) routine.
"""
function gmst(::IAU1980Model, t::Number)

    # Coefficients for the IAU 1982 GMST-UT1 model. The first component has been adjusted 
    # of 12 hours because UT1 starts at noon.
    A = -19089.45159
    B = 8640184.812866
    C = 0.093104
    D = -6.2e-6

    # Fractional part of UT1, in seconds
    f = Tempo.CENTURY2SEC*mod(t, 1)

    # seconds to radians conversion 
    s2r = 7.272205216643039903848712e-5
    
    # Compute GMST
    return mod2pi(s2r*(@evalpoly(t, A, B, C, D) + f))

end


"""
Compute the Greenwich Apparent Sidereal Time (GAST), in radians, given time `t` as `TT` 
Julian centuries since `J2000`
"""
function gast(m::IAU2006Model, t::Number, θ::Number)

    # TODO: or maybe its the transpose?
    NPB = orient_bias_precession_nutation(m, t)
    
    # Get the CIO locator 
    s = cip_coords(m, t)

    # Compute the equations of the origins
    eors = origins_equation()  

    return θ - eors
end


function gast(m::IAU2000Model, t::Number, θ::Number)

    # Compute the Greenwich Mean Sidereal time 
    gmst = gmst(m, θ, t)

    # Compute the equations of the equinoxes
    ee = equinoxes_equation(m, t)

    # Comput GAST, in radians
    return mod2pi(gmst + ee)

end


function gast(::IAU1980Model)

    # Compute the Greenwich Mean Sidereal time 
    gmst = gmst(m, t)

    # Compute the equations of the equinoxes
    ee = equinoxes_equation(m, t)

    # Comput GAST, in radians
    return mod2pi(gmst + ee)

end


# Time expressed in TDB julian centuries since J2000
function equinoxes_equation(m::IAU2006Model, t::Number)

    # Compute GMST and GAST, in radians 
    gmst = gmst(m, t, θ)
    gast = gast(m, t, θ)

    # Equation of the equinoxes
    return mod2pi(gast - gmst)

end

# Time expressed in TT julian centuries since J2000
function equinoxes_equation(m::IAU2000Model, t::Number)

    # Compute precession-rate adjustments 

    # Compute the mean obliquity 
    ϵₐ = orient_obliquity(m, t)

    # Nutation in logitutude 

    # Equation of the equinoxes
    # TODO: here you have to implemetn a damned series 
    return Δψ*cos(ϵₐ)

end

# Time expressed in TDB julian centuries since J2000
function equinoxes_equation(m::IAU1980Model, t::Number)

    # Longitude of the mean ascending node of the lunar orbit on the ecliptic, 
    # measured from the mean equinox of date, in radians
    ω = mod2pi(arcsec2rad(@evalpoly(t, -482890.539, 7.455, 0.008)) + 2π*mod(-5t, 1))

    # Compute the nutation components in longitude and obliquity
    # TODO: missing orient_nutation per il 1980 missa 

    # Compute the mean obliquity 
    ϵₐ = orient_obliquity(m, t)

    # Equations of the equinoxes 
    return Δψ*cos(ϵₐ) + arcsec2rad(0.00264*sin(ω) + 0.000063*sin(2ω))

end


function origins_equation(::IAU2006Model, t::Number)

end