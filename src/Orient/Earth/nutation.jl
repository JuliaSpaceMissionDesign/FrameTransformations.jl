export orient_nutation

include("constants/nut2000a.jl");
include("constants/nut2000b.jl");

# Arcseconds in a full circle
const ARCSECTURN = 1296000.0

nutation00(::IAU2006Model, ::Number, ::FundamentalArguments) = () 
build_nutation_series(:nutation00, :IAU2006A,  NUTATION_2000Aψ, NUTATION_2000Aϵ)
build_nutation_series(:nutation00, :IAU2006B,  NUTATION_2000Bψ, NUTATION_2000Bϵ)


"""
    orient_nutation(m::IAU2006Model, t::Number)

Compute the longitude and obliquity components of the IAU 2006/2000 A/B
precession-nutation model, in radians.

### Inputs 
- `m` -- IAU 2006 Precession-Nutation Model 
- `t` -- Terrestrial Time `TT` in Julian Centuries since J2000. 

### Notes 
- Due to their theoretical basis, the original developments required `t` expressed 
as TDB. However, in practice, it is usually more convenient to use Terrestrial Time (TT)
as it makes no significant differences (< 0.01 μas) in the final result.

- For the IAU 2006A model, the function strictly follows the SOFA implementation. 
It first computes the IAU 2000A nutation, then applies adjustments for the 
consequences of the change in obliquity from the IAU 1980 ecliptic to the IAU 2006 
ecliptic and (ii) for the secular variation in the Earth's dynamical form factor J2. 
These corrections ensure that the IAU 2000A nutation is consistent with the IAU 2006 
precession model. Please note that the coefficients available on the IERS tables 
already include those corrections, and are retrieved by multiplied the amplitudes 
of the SOFA nutation in longitude coefficients by 1.00000047. 

- The computation of the free-core nutation and time dependent effects are excluded 
from this model. To achieve the < 1μas accuracy with the IAU 2006A precession-nutation model, 
such effects must be included a-posteriori (through dX and dY) using the IERS EOP data.

- For the IAU 2000B model, the nutation series is truncated from nearly 1400 terms 
to only 77, yet it still delivers results of 1 mas accuracy at present epochs. In 
particular, it delivers a pole accurate to 1 mas from 1900 to 2100 (only very 
occasionally just outside 1 mas). The coefficients are taken from SOFA's implementation, 
which slighlty differ from those reported in McCarthy and Luzum (2003). Comparisons with 
IAU 2006A show that the SOFA version between 1995 and 2050 delivers 0.283 mas RMSE 
(0.994 mas in the worst case), whereas the IERS Conventions website version delivers 
0.312 RMSE (1.125 mas in the worst case).

- The IAU 2000B model includes planetary bias terms that compensate for long-period 
nutations. These amplitudes used in this implementation are optimised for a _rigorous_
method, where frame bias, precession and nutation are applied separately and in that order 
(see SOFA's documentation for further insights).

- A simplified version of the Fundamental Arguments, taken from Simon et al (1994) 
is exploited for IAU2006B as the erorr introduced is below the model accuracy ( ~0.1 mas).

### References 
- Luzum, B. and Petit G. (2012), _The IERS Conventions (2010)_, 
[IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- Wallace P. T. and Capitaine N. (2006), _Precession-nutation procedures consistent with 
IAU 2006 resolutions_, [DOI: 10.1051/0004-6361:20065897](https://www.aanda.org/articles/aa/abs/2006/45/aa5897-06/aa5897-06.html) 
- Simon, J. et al., (1994), _Numerical expressions for precession formulae and mean elements for the 
Moon and the planets_.
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/pfw06.c) library
"""
function orient_nutation(m::IAU2006A, t::Number)

    # Computes Fundamental Arguments from IERS 2003
    fa = FundamentalArguments(t)

    # Computes IAU 2000A nutation components from luni-solar 
    # and planetary terms of the Mathews et al. (2002) series
    δψₐ, δϵₐ = nutation00(m, t, fa)

    # Factor correcting the secular variation of J2 
    f = -2.7774e-6*t # t = TT 

    # Applies P03 Nutation Corrections from WC06 (2006)
    Δψ = (1 + f + 0.4697e-6)*δψₐ
    Δϵ = (1 + f)*δϵₐ

    return Δψ, Δϵ
end

function orient_nutation(m::IAU2006B, t::Number)

    # Computes only Luni-Solar Fundamental Arguments 
    fa = FundamentalArguments(LuniSolarArguments(t, m)..., 0., 0., 
                              0., 0., 0., 0., 0., 0., 0.)

    # Computes luni-solar nutation contributions 
    δψ_ls, δϵ_ls = nutation00(m, t, fa)

    # Adds offset to account for truncated planetary contributions 
    δψ_pl = arcsec2rad(-0.135 * 1e-3)
    δϵ_pl = arcsec2rad( 0.388 * 1e-3)

    return δψ_ls + δψ_pl, δϵ_ls + δϵ_pl
end

