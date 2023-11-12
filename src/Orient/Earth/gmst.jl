export gmst

"""
    gmst(θ::Number, t::Number)

Compute the Greenwich Mean Sidereal Time from the Earth rotation angle `Θ` and the epoch  `t` 
expressed in `TT` Julian centuries since `J2000`. 
"""
function gmst(θ::Number, t::Number)
    c₀ = 0.014506
    c₁ = 4612.156534
    c₂ = 1.3915817
    c₃ = -0.00000044
    c₄ = -0.000029956
    c₅ = -0.0000000368
    return 86400 * θ + @evalpoly(t, c₀, c₁, c₂, c₃, c₄, c₅)/15
end