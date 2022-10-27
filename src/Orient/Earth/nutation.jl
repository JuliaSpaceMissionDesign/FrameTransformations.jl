export orient_nutation

struct NutationLuniSolar
    # Coefficients of l,l',F,D,Om
    nl::Int
    nlp::Int
    nf::Int
    nd::Int
    nom::Int
    # longitude sin, t*sin, cos coefficients
    sp::Float64
    spt::Float64
    cp::Float64
    # obliquity cos, t*cos, sin coefficients
    ce::Float64
    cet::Float64
    se::Float64
end

struct NutationPlanetary
    # Coefficients of l, F, D and Omega
    nl::Int
    nf::Int
    nd::Int
    nom::Int
    # Coefficients of planetary longitudes
    nme::Int
    nve::Int
    nea::Int
    nma::Int
    nju::Int
    nsa::Int
    nur::Int
    nne::Int
    # Coefficient of general precession
    npa::Int
    # Longitude sin, cos coefficients
    sp::Int
    cp::Int
    # Obliquity sin, cos coefficients
    se::Int
    ce::Int
end

include(joinpath("constants", "nutation2000b.jl"))

# Arcseconds in a full circle
const ARCSECTURN = 1296000.0


"""
    nutation(::IAU2006B, t::N) where {N<:Number}

Nutation, IAU 2000B model.

### Inputs 
- `t`-- `TT` centuries since J2000

### Output 
`dψ, dε` nutation, lunisolar + planetary terms.

### Notes
The nutation components in longitude and obliquity are in radians and with 
respect to the equinox and ecliptic of date. The obliquity at J2000.0 is assumed 
to be the Lieske et al. (1977) value of 84381.448 arcsec.

The nutation model consists only of luni-solar terms, but includes also a fixed 
offset which compensates for certain long-period planetary terms.

The IAU 2000B model (McCarthy & Luzum 2003) contains only 77 terms, plus additional 
simplifications, yet still deliver results of **1 mas accuracy at present epochs**.

### References 
- Lieske, J.H., Lederle, T., Fricke, W., Morando, B., "Expressions for the 
    precession quantities based upon the IAU /1976/ system of astronomical 
    constants", Astron.Astrophys. 58, 1-2, 1-16. (1977)
- McCarthy, D.D. & Luzum, B.J., "An abridged model of the precession-nutation 
    of the celestial pole", Cel.Mech.Dyn.Astron. 85, 37-49 (2003)
- Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M., Francou, G., 
    Laskar, J., Astron.Astrophys. 282, 663-683 (1994)
- [AstroBase](https://github.com/JuliaAstro/AstroBase.jl/blob/master/src/EarthAttitude/nutation.jl) software library
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/nut00b.c) software library
"""
function orient_nutation(::IAU2006B, t::N) where {N<:Number} 
    # Fundamental (Delaunay) arguments from Simon et al. (1994)
    # Mean anomaly of the Moon
    el = mod(485868.249036+(1717915923.2178) * t, ARCSECTURN) |> arcsec2rad

    # Mean anomaly of the Sun. 
    elp = mod(1287104.79305+(129596581.0481) * t, ARCSECTURN) |> arcsec2rad

    # Mean argument of the latitude of the Moon. 
    f = mod(335779.526232+(1739527262.8478) * t, ARCSECTURN) |> arcsec2rad

    # Mean elongation of the Moon from the Sun. 
    d = mod(1072260.70369+ (1602961601.2090) * t, ARCSECTURN) |> arcsec2rad

    # Mean longitude of the ascending node of the Moon. 
    om = mod(450160.398036 + (-6962890.5431) * t, ARCSECTURN) |> arcsec2rad

    # Initialize the nutation values.
    dp = 0.0 
    dq = 0.0

    for x in Iterators.reverse(NUTATION_2000B)
        arg = mod2pi(x.nl * el + x.nlp * elp + x.nf * f + x.nd * d + x.nom * om)
        sarg, carg = sincos(arg)

        dp += (x.sp + x.spt * t) * sarg + x.cp * carg
        dq += (x.ce + x.cet * t) * carg + x.se * sarg
    end

    # Luni-solar terms
    δψ_ls = arcsec2rad(dp * 1e-7)
    δϵ_ls = arcsec2rad(dq * 1e-7)

    # Fixed offset to correct for missing terms in truncated series.
    δψ_pl = arcsec2rad(-0.135 * 1e-3)
    δϵ_pl = arcsec2rad(0.388 * 1e-3)

    return δψ_ls + δψ_pl, δϵ_ls + δϵ_pl
end