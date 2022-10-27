const EQ_PREC_POL = (
    (5453.282155, 0.4252841, -0.00037173, -0.000000152),
    (-73750.930350, -0.7675452, -0.00018725, 0.000000231)
)

const EQ_PREC_PER = (
    ( 256.75, -819.940624,75004.344875,81491.287984, 1558.515853),
    ( 708.15,-8444.676815,  624.033993,  787.163481, 7774.939698),
    ( 274.20, 2600.009459, 1251.136893, 1251.296102,-2219.534038),
    ( 241.45, 2755.175630,-1102.212834,-1257.950837,-2523.969396),
    (2309.00, -167.659835,-2660.664980,-2966.799730,  247.850422),
    ( 492.20,  871.855056,  699.291817,  639.744522, -846.485643),
    ( 396.10,   44.769698,  153.167220,  131.600209,-1393.124055),
    ( 288.90, -512.313065, -950.865637, -445.040117,  368.526116),
    ( 231.10, -819.415595,  499.754645,  584.522874,  749.045012),
    (1610.00, -538.071099, -145.188210,  -89.756563,  444.704518),
    ( 620.00, -189.793622,  558.116553,  524.429630,  235.934465),
    ( 157.87, -402.922932,  -23.923029,  -13.549067,  374.049623),
    ( 220.30,  179.516345, -165.405086, -210.157124, -171.330180),
    (1200.00,   -9.814756,    9.344131,  -44.919798,  -22.899655)
)


"""
    lt_equator_precession(j2000ttc::N) where {N<:Number}

Long-term precession of the equator.

!!! note 
    The returned vector is with respect to the J2000.0 mean equator and equinox.

### Inputs 
- `j2000ttc` -- `TT` centuries since J2000

### Output 
Equator pole unit vector.

The Vondrak et al. (2011, 2012) 400 millennia precession model agrees with the 
IAU 2006 precession at J2000.0 and stays within 100 microarcseconds during the 
20th and 21st centuries.  It is accurate to a few arcseconds throughout the 
historical period, worsening to a few tenths of a degree at the end of the  
+/- 200,000 year time span.

### References 
- Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession expressions, 
    valid for long time intervals, Astron.Astrophys. 534, A22
- Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession expressions, 
    valid for long time intervals (Corrigendum), Astron.Astrophys. 541, C1
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/ltpequ.c) library
"""
function lt_equator_precession(j2000ttc::N) where {N<:Number}
    # Initialize X and Y accumulators
    x = 0.0
    y = 0.0

    # Periodic terms 
    w = 2π * j2000ttc
    for i in eachindex(EQ_PREC_PER)
        a = w/EQ_PREC_PER[i][1]
        s, c = sincos(a)
        x += c*EQ_PREC_PER[i][2] + s*EQ_PREC_PER[i][4]
        y += c*EQ_PREC_PER[i][3] + s*EQ_PREC_PER[i][5]
    end

    # Polynomial terms
    w = 1.0
    for i in 1:length(EQ_PREC_POL[1])
        x += EQ_PREC_POL[1][i]*w
        y += EQ_PREC_POL[2][i]*w
        w *= j2000ttc
    end

    x = arcsec2rad(x)
    y = arcsec2rad(y)
    w = 1.0 - x*x - y*y 

    return SA[
        x, 
        y, 
        w < 0.0 ? 0.0 : sqrt(w)
    ]
end

const EC_PREC_POL = (
    (5851.607687, -0.1189000, -0.00028913, 0.000000101),
    (-1600.886300, 1.1689818, -0.00000020, -0.000000437)
)

const EC_PREC_PER = (
    ( 708.15,-5486.751211,-684.661560,  667.666730,-5523.863691),
    (2309.00,  -17.127623,2446.283880,-2354.886252, -549.747450),
    (1620.00, -617.517403, 399.671049, -428.152441, -310.998056),
    ( 492.20,  413.442940,-356.652376,  376.202861,  421.535876),
    (1183.00,   78.614193,-186.387003,  184.778874,  -36.776172),
    ( 622.00, -180.732815,-316.800070,  335.321713, -145.278396),
    ( 882.00,  -87.676083, 198.296701, -185.138669,  -34.744450),
    ( 547.00,   46.140315, 101.135679, -120.972830,   22.885731)
)

"""
    lt_ecliptic_precession(j2000tt::N) where {N<:Number}

Long-term ecliptic of the equator.

### Inputs 
- `j2000ttc` -- `TT` centuries since J2000.

### Output 
Ecliptic pole unit vector.

!!! note 
    The returned vector is with respect to the J2000.0 mean equator and equinox.

The Vondrak et al. (2011, 2012) 400 millennia precession model agrees with the 
IAU 2006 precession at J2000.0 and stays within 100 microarcseconds during the 
20th and 21st centuries.  It is accurate to a few arcseconds throughout the 
historical period, worsening to a few tenths of a degree at the end of the  
+/- 200,000 year time span.

### References 
- Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession expressions, 
    valid for long time intervals, Astron.Astrophys. 534, A22
- Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession expressions, 
    valid for long time intervals (Corrigendum), Astron.Astrophys. 541, C1
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/ltpequ.c) library
"""
function lt_ecliptic_precession(j2000ttc::N) where {N<:Number}
    # Initialize X and Y accumulators
    p = 0.0
    q = 0.0

    # Periodic terms 
    w = 2π * j2000ttc
    for i in eachindex(EC_PREC_PER)
        a = w/EC_PREC_PER[i][1]
        s, c = sincos(a)
        p += c*EC_PREC_PER[i][2] + s*EC_PREC_PER[i][4]
        q += c*EC_PREC_PER[i][3] + s*EC_PREC_PER[i][5]
    end

    # Polynomial terms
    w = 1.0
    for i in 1:length(EC_PREC_POL[1])
        p += EC_PREC_POL[1][i]*w
        q += EC_PREC_POL[2][i]*w
        w *= j2000ttc
    end

    # Form the ecliptic pole vector
    p = arcsec2rad(p)
    q = arcsec2rad(q)
    w = 1.0 - p*p + q*q
    ϵ0 = 0.4090926006005829 # arcsec2rad(84381.406)
    s, c = sincos(ϵ0)
    
    return SA[
        p, 
        - q*c - w*s,
        - q*s + w*c
    ]
end

"""
    lt_precession(j2000ttc::N) where {N<:Number}

Long term precession matrix.

### Inputs 
- `j2000ttc` -- `TT` centuries since J2000.

### Output
Precession matrix. 

The matrix is in the sense P(date) = RP x P(J2000), where P(J2000) is a vector 
with respect to the J2000.0 mean equator and equinox and P(date) is the same 
vector with respect to the equator and equinox of the date.

The Vondrak et al. (2011, 2012) 400 millennia precession model agrees with the 
IAU 2006 precession at J2000.0 and stays within 100 microarcseconds during the 
20th and 21st centuries.  It is accurate to a few arcseconds throughout the 
historical period, worsening to a few tenths of a degree at the end of the  
+/- 200,000 year time span.

### References 
- Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession expressions, 
    valid for long time intervals, Astron.Astrophys. 534, A22
- Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession expressions, 
    valid for long time intervals (Corrigendum), Astron.Astrophys. 541, C1
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/ltp.c) library
"""
function lt_precession(j2000ttc::N) where {N<:Number}
    # Equator pole (bottom row of matrix)
    peqr = lt_equator_precession(j2000ttc)

    #  Ecliptic pole
    pecl = lt_ecliptic_precession(j2000ttc)

    # Equinox (top row of matrix).
    v = cross(peqr, pecl)
    eqx = v./norm(v)

    # Middle row of matrix
    v = cross(peqr, eqx)

    return SMatrix{3, 3}(
        eqx[1],     eqx[2],     eqx[3],
        v[1],         v[2],       v[3],
        peqr[1],   peqr[2],    peqr[3]
    )'
end

"""
    lt_precession_bias(j2000ttc::N) where {N<:Number}

Long term precession matrix, including ICRS frame bias.

### Inputs 
- `j2000ttc` -- `TT` centuries since J2000.

### Output
Precession matrix. 

The matrix is in the sense P(date) = RP x P(ICRS), where P(ICRS) is a vector in 
the Geocentric Celestial Reference and P(date) is the vector with respect 
to the Celestial Intermediate Reference System at that date but with nutation 
neglected.

!!! note 
    A first order frame bias formulation is used, of sub-microarcsecond 
    accuracy compared with a full 3D rotation.

The Vondrak et al. (2011, 2012) 400 millennia precession model agrees with the 
IAU 2006 precession at J2000.0 and stays within 100 microarcseconds during the 
20th and 21st centuries.  It is accurate to a few arcseconds throughout the 
historical period, worsening to a few tenths of a degree at the end of the  
+/- 200,000 year time span.

### References 
- Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession expressions, 
    valid for long time intervals, Astron.Astrophys. 534, A22
- Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession expressions, 
    valid for long time intervals (Corrigendum), Astron.Astrophys. 541, C1
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/ltp.c) library
"""
function lt_precession_bias(j2000ttc::N) where {N<:Number}
    # Precession matrix
    P = lt_precession(j2000ttc)

    # Frame bias (IERS Conventions 2010, Eqs. 5.21 and 5.33)
    dx = -8.056148938997159e-8 # arcsec2rad(-0.016617)
    de = -3.3060414542221477e-8 # arcsec2rad(-0.0068192)
    dr = -7.078279744199226e-8 # arcsec2rad(-0.0146)

    return SMatrix{3, 3}(
        P[1,1]-P[1,2]*dr+P[1,3]*dx, P[2,1]-P[2,2]*dr+P[2,3]*dx, P[3,1]-P[3,2]*dr+P[3,3]*dx,
        P[1,1]*dr+P[1,2]+P[1,3]*de, P[2,1]*dr+P[2,2]+P[2,3]*de, P[3,1]*dr+P[3,2]+P[3,3]*de,
        P[1,1]*dx-P[1,2]*de+P[1,3], P[2,1]*dx-P[2,2]*de+P[2,3], P[3,1]*dx-P[3,2]*de+P[3,3]
    )
end