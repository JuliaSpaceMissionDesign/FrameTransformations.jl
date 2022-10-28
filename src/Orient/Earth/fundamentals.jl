
struct FundamentalArguments{N <: Number}

    # Planetary arguments (mean longitutes)
    λ_Me::N 
    λ_Ve::N 
    λ_Ea::N 
    λ_Ma::N 
    λ_Ju::N 
    λ_Sa::N 
    λ_Ur::N
    λ_Ne::N

    # Moon-Sun Delaunay Arguments
    uₘ::N # moon latitude
    Dₛ::N # moon elongation
    Ωₘ::N # moon raan 
    Mₐ::N # moon anomaly 
    Sₐ::N # sun anomaly

    pₐ::N # general accumulated precession

end

# Computes Fundamental Arguments at epoch t
function FundamentalArguments(t::Number)
    FundamentalArguments( 
        fa_mlon_mercury(t), 
        fa_mlon_venus(t), 
        fa_mlon_earth(t), 
        fa_mlon_mars(t),
        fa_mlon_jupiter(t),
        fa_mlon_saturn(t), 
        fa_mlon_uranus(t), 
        fa_mlon_neptune(t), 
        fa_mlat_moon(t),
        fa_melo_moon(t),
        fa_mlon_moon(t), 
        fa_mano_moon(t),
        fa_mano_sun(t),
        fa_precession(t)
    )
end

# Computes mean anomaly of the moon 
function fa_mano_moon(t::Number) 
    mod2pi(@evalpoly(t, 485868.249036, 
                        1717915923.2178, 
                        31.8792, 
                        0.051635, 
                        -0.00024470) |> arcsec2rad)
end


# Computes mean anomaly of the sun
function fa_mano_sun(t::Number) 
    mod2pi(@evalpoly(t, 1287104.793048, 
                        129596581.0481, 
                        -0.5532,
                        0.000136, 
                        -0.00001149) |> arcsec2rad)
end

# Computes mean longitude of the moon's ascending node  
function fa_mlon_moon(t::Number) 
    mod2pi(@evalpoly(t, 450160.398036, 
                        -6962890.54311, 
                        7.4722,
                        0.007702, 
                        -0.00005939) |> arcsec2rad)
end

# Computes mean argument of latitude of the Moon
function fa_mlat_moon(t::Number)
    mod2pi(@evalpoly(t, 335779.526232, 
                        1739527262.8478, 
                        -12.7512,
                        -0.001037, 
                        -0.00000417) |> arcsec2rad)
end

# Computes mean elongation of the Moon from the Sun 
function fa_melo_moon(t::Number)
    mod2pi(@evalpoly(t, 1072260.70369 , 
                        1602961601.2090, 
                        -6.3706,
                        -0.006593, 
                        -0.00003169) |> arcsec2rad)
end

# Compute general accumulated precession in longitude 
function fa_precession(t::Number)
    mod2pi(@evalpoly(t, 0, 0.024381750, 0.00000538691))
end

# Compute mean heliocentric longitude of Mercury 
function fa_mlon_mercury(t::Number)
    mod2pi(@evalpoly(t, 4.402608842, 2608.7903141574))
end

# Compute mean heliocentric longitude of Venus 
function fa_mlon_venus(t::Number)
    mod2pi(@evalpoly(t, 3.176146697, 1021.3285546211))
end

# Compute mean heliocentric longitude of Earth 
function fa_mlon_earth(t::Number)
    mod2pi(@evalpoly(t, 1.753470314, 628.3075849991))
end

# Compute mean heliocentric longitude of Mars 
function fa_mlon_mars(t::Number)
    mod2pi(@evalpoly(t, 6.203480913, 334.0612426700))
end

# Compute mean heliocentric longitude of Jupiter 
function fa_mlon_jupiter(t::Number)
    mod2pi(@evalpoly(t, 0.599546497, 52.9690962641))
end

# Compute mean heliocentric longitude of Saturn 
function fa_mlon_saturn(t::Number)
    mod2pi(@evalpoly(t, 0.874016757, 21.3299104960))
end

# Compute mean heliocentric longitude of Uranus 
function fa_mlon_uranus(t::Number)
    mod2pi(@evalpoly(t, 5.481293872, 7.4781598567))
end

# Compute mean heliocentric longitude of Neptune 
function fa_mlon_neptune(t::Number)
    mod2pi(@evalpoly(t, 5.311886287, 3.8133035638))
end
