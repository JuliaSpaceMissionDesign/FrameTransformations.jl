export orient_nutation

include("constants/nut2000a.jl");
include("constants/nut2000b.jl");

# Arcseconds in a full circle
const ARCSECTURN = 1296000.0

nutation00(::IAU2006Model, ::Number, ::FundamentalArguments) = () 
build_nutation_series(:nutation00, :IAU2006A,  NUTATION_2000Aψ, NUTATION_2000Aϵ)
build_nutation_series(:nutation00, :IAU2006B,  NUTATION_2000Bψ, NUTATION_2000Bϵ)

# ### Notes
# This implementation follows SOFA routines, where even though  `t` is strictly 
# expressed as TDB, it is usually more convenient to use Terrestrial Time (TT), 
# which makes no significant differences (< 0.01 μas) in the final result.

function orient_nutation(m::IAU2006A, t::Number, fa::FundamentalArguments)

    # t should be in Terrestrial Time (TT)

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

function orient_nutation(m::IAU2006A, t::Number)
    fa = FundamentalArguments(t)
    orient_nutation(m, t, fa)
end


function orient_nutation(m::IAU2006B, t::Number, fa::FundamentalArguments)

    # Computes luni-solar nutation contributions 
    δψ_ls, δϵ_ls = nutation00(m, t, fa)

    # Adds offset to account for truncated planetary contributions 
    δψ_pl = arcsec2rad(-0.135 * 1e-3)
    δϵ_pl = arcsec2rad( 0.388 * 1e-3)

    return δψ_ls + δψ_pl, δϵ_ls + δϵ_pl
end


function orient_nutation(m::IAU2006B, t::Number)

    # questa utilizza una versione semplificata degli FA, tanto l'errore introdotto è 
    # minore della precisione del modello, sui circa 0.1 mas

    # Mean anomalies of the Moon and Sun
    Mₐ = mod(485868.249036 + 1717915923.2178t, ARCSECTURN) |> arcsec2rad
    Sₐ = mod(1287104.79305 + 129596581.0481t, ARCSECTURN) |> arcsec2rad

    # Mean argument of the latitude of the Moon. 
    uₘ = mod(335779.526232 + 1739527262.8478t, ARCSECTURN) |> arcsec2rad

    # Mean elongation of the Moon from the Sun. 
    Dₛ = mod(1072260.70369 + 1602961601.2090t, ARCSECTURN) |> arcsec2rad

    # Mean longitude of the ascending node of the Moon. 
    Ωₘ = mod(450160.398036 -6962890.5431t, ARCSECTURN) |> arcsec2rad

    fa = FundamentalArguments(0., 0., 0., 0., 0., 0., 0., 
                              0., uₘ, Dₛ, Ωₘ, Mₐ, Sₐ, 0.)

    orient_nutation(m, t, fa)
end

