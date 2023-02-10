
"""
	FundamentalArguments{N <: Number}

Type storing the fundamental luni-solar and planetary arguments.

### Fields 
- `Mₐ` -- Mean anomaly of the Moon 
- `Sₐ` -- Mean anomaly of the Sun
- `uₘ` -- Mean longitude of the Moon minus mean longitude of the ascending node `F`
- `Dₛ` -- Mean elongation of the Moon from the Sun 
- `Ωₘ` -- Mean longitude of the Moon's ascending node
- `λ_Me` -- Mercury's mean heliocentric longitude
- `λ_Ve` -- Venus's mean heliocentric longitude
- `λ_Ea` -- Earth's mean heliocentric longitude
- `λ_Ma` -- Mars's mean heliocentric longitude
- `λ_Ju` -- Jupiter's mean heliocentric longitude
- `λ_Sa` -- Saturn's mean heliocentric longitude
- `λ_Ur` -- Uranus's mean heliocentric longitude
- `λ_Ne` -- Neptune's mean heliocentric longitude
- `pₐ` -- General precession in longitude
"""
struct FundamentalArguments{N <: Number}

	# Luni-Solar Arguments
	Mₐ::N  
	Sₐ::N 
	uₘ::N 
	Dₛ::N 
	Ωₘ::N  

	# Planetary Mean Heliocentric Longitudes
	λ_Me::N 
	λ_Ve::N 
	λ_Ea::N 
	λ_Ma::N 
	λ_Ju::N 
	λ_Sa::N 
	λ_Ur::N
	λ_Ne::N

	pₐ::N # General Accumulated Precession

end

"""
	FundamentalArguments(t::Number, m::IAUModel=iau2006a)

Compute the fundamental luni-solar and planetary arguments, in radians, associated to the 
given IAU model `m`, at time `t` expressed in `TDB` Julian centuries since [`J2000`](@ref).

The available IAU model options are: [`iau2006a`](@ref) and [`iau2000b`](@ref).

!!! note 
	Only the planetary arguments are affected by the IAU Model choice.
"""
function FundamentalArguments(t::Number, m::IAUModel=iau2006a)
	FundamentalArguments(LuniSolarArguments(t, m)..., PlanetaryArguments(t)...)
end

"""
	LuniSolarArguments(t::Number, m::IAUModel)

Compute the fundamental (Delaunay) Luni-Solar arguments, in radians, associated to the 
desired IAU model `m`, at time `t` expressed in `TDB` Julian centuries since [`J2000`](@ref).

The returned values depend on the input model as follows: 

- **IAU2006A**: the Delaunay expressions are taken from the IERS 2010 Conventions.

- **IAU2000B**: the expressions are taken from Simon et al. (1994), following ERFA's
    implementation of nut00b.c


# Outputs 
- `Mₐ` -- Mean anomaly of the Moon
- `Sₐ` -- Mean anomaly of the Sun
- `uₘ` -- Mean longitude of the Moon minus mean longitude of the
          ascending node `F`
- `Dₛ` -- Mean elongation of the Moon from the Sun
- `Ωₘ` -- Mean longitude of the Moon's ascending node

### References 
- Luzum, B. and Petit G. (2012). The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- Simon, J. et al., (1994), Numerical expressions for precession formulae and mean elements for the 
  Moon and the planets.
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/nut00b.c) software library
"""
function LuniSolarArguments(t::Number, ::Union{<:IAU2000A, <:IAU2006A}) 
	Mₐ = fa_mano_moon(t)
	Sₐ = fa_mano_sun(t)
	uₘ = fa_mlat_moon(t) 
	Dₛ = fa_melo_moon(t)
	Ωₘ = fa_mlon_moon(t)

	return Mₐ, Sₐ, uₘ, Dₛ, Ωₘ
end

function LuniSolarArguments(t::Number, ::Union{<:IAU2000B, <:IAU2006B, <:CPNModel})

	# Approximated values as linear functions of time consistent with IAU2000B
	# Mean anomalies of the Moon and Sun
	Mₐ = mod(485868.249036 + 1717915923.2178t, ARCSECTURN) |> arcsec2rad
	Sₐ = mod(1287104.79305 + 129596581.0481t, ARCSECTURN) |> arcsec2rad

	# Mean argument of the latitude of the Moon. 
	uₘ = mod(335779.526232 + 1739527262.8478t, ARCSECTURN) |> arcsec2rad

	# Mean elongation of the Moon from the Sun. 
	Dₛ = mod(1072260.70369 + 1602961601.2090t, ARCSECTURN) |> arcsec2rad

	# Mean longitude of the ascending node of the Moon. 
	Ωₘ = mod(450160.398036 -6962890.5431t, ARCSECTURN) |> arcsec2rad

	return Mₐ, Sₐ, uₘ, Dₛ, Ωₘ
end


"""
	PlanetaryArguments(t::Number)

Compute the fundamental planetary arguments and the general precession, in radians, 
at time `t` expressed in `TDB` Julian centuries since [`J2000`](@ref).

### Outputs 
- `λ☿` -- Mercury's mean heliocentric longitude.
- `λ♀` -- Venus's mean heliocentric longitude.
- `λe` -- Earth's mean heliocentric longitude.
- `λ♂` -- Mars's mean heliocentric longitude.
- `λ♃` -- Jupiter's mean heliocentric longitude.
- `λ♄` -- Saturn's mean heliocentric longitude.
- `λ⛢` -- Uranus's mean heliocentric longitude.
- `λ♆` -- Neptune's mean heliocentric longitude.
- `pλ` -- General precession in longitude.

### References 
- Luzum, B. and Petit G. (2012). The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
"""
function PlanetaryArguments(t::Number)
	λ_Me = fa_mlon_mercury(t)
	λ_Ve = fa_mlon_venus(t)
	λ_Ea = fa_mlon_earth(t)
	λ_Ma = fa_mlon_mars(t)
	λ_Ju = fa_mlon_jupiter(t)
	λ_Sa = fa_mlon_saturn(t)
	λ_Ur = fa_mlon_uranus(t)
	λ_Ne = fa_mlon_neptune(t)
	pₐ   = fa_precession(t)

	return λ_Me, λ_Ve, λ_Ea, λ_Ma, λ_Ju, λ_Sa, λ_Ur, λ_Ne, pₐ
end


"""
	fa_mano_moon(t::Number) 

Return the mean anomaly of the Moon `Mₐ`, in radians, at time `t` expressed in `TDB` Julian 
centuries since [`J2000`](@ref).

### References 
- Luzum, B. and Petit G. (2012). The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/fal03.c) software library
"""
function fa_mano_moon(t::Number) 
	mod2pi(@evalpoly(t, 485868.249036, 
						1717915923.2178, 
						31.8792, 
						0.051635, 
						-0.00024470)*π/648000)
end


"""
	fa_mano_sun(t::Number) 

Return the mean anomaly of the Sun `Sₐ` in radians, at time `t` expressed in `TDB` Julian 
centuries since [`J2000`](@ref).

### References 
- Luzum, B. and Petit G. (2012). The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/falp03.c) software library
"""
function fa_mano_sun(t::Number) 
	mod2pi(@evalpoly(t, 1287104.793048, 
						129596581.0481, 
						-0.5532,
						0.000136, 
						-0.00001149)*π/648000)
end


"""
	fa_mlon_moon(t::Number) 

Return the mean longitude of the Moon's ascending node `Ω` in radians, at time `t` expressed
in `TDB` Julian centuries since [`J2000`](@ref).

### References 
- Luzum, B. and Petit G. (2012). The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/faom03.c) software library
"""
function fa_mlon_moon(t::Number) 
	mod2pi(@evalpoly(t, 450160.398036, 
						-6962890.5431, 
						7.4722,
						0.007702, 
						-0.00005939)*π/648000)
end


"""
	fa_mlat_moon(t::Number) 

Return the mean longitude of the Moon minus the mean longitude of the ascending node `F` in
radians, at time `t` expressed in `TDB` Julian centuries since [`J2000`](@ref).  

### References 
- Luzum, B. and Petit G. (2012). The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/faf03.c) software library
"""
function fa_mlat_moon(t::Number)
	mod2pi(@evalpoly(t, 335779.526232, 
						1739527262.8478, 
						-12.7512,
						-0.001037, 
						+0.00000417)*π/648000)
end


"""
	fa_melo_moon(t::Number) 

Return the mean elongation of the Moon from the Sun `D` in radians, at time `t` expressed in 
`TDB` Julian centuries since [`J2000`](@ref).

### References 
- Luzum, B. and Petit G. (2012). The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/fad03.c) software library
"""
function fa_melo_moon(t::Number)
	mod2pi(@evalpoly(t, 1072260.703692 , 
						1602961601.2090, 
						-6.3706,
						+0.006593, 
						-0.00003169)*π/648000)
end


"""
  	fa_precession(t::Number) 

Return the general accumulated precession in longitude `pₐ` in radians, at time `t` 
expressed in `TDB` Julian centuries since [`J2000`](@ref). 

### References 
- Luzum, B. and Petit G. (2012). The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/fapa03.c) software library
"""
function fa_precession(t::Number)
  	mod2pi(@evalpoly(t, 0, 0.024381750, 0.00000538691))
end


"""
	fa_mlon_mercury(t::Number) 

Return the mean heliocentric longitude of Mercury in radians, at time `t` expressed in `TDB` 
Julian centuries since [`J2000`](@ref). 

### References 
- Luzum, B. and Petit G. (2012). The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/fame03.c) software library
"""
function fa_mlon_mercury(t::Number)
  	mod2pi(@evalpoly(t, 4.402608842, 2608.7903141574))
end


"""
	fa_mlon_venus(t::Number) 

Return the mean heliocentric longitude of Venus in radians, at time `t` expressed in `TDB` 
Julian centuries since [`J2000`](@ref). 

### References 
- Luzum, B. and Petit G. (2012). The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/fave03.c) software library
"""
function fa_mlon_venus(t::Number)
  	mod2pi(@evalpoly(t, 3.176146697, 1021.3285546211))
end


"""
	fa_mlon_earth(t::Number) 

Return the mean heliocentric longitude of Earth in radians, at time `t` expressed in `TDB` 
Julian centuries since [`J2000`](@ref). 

### References 
- Luzum, B. and Petit G. (2012). The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/fae03.c) software library
"""
function fa_mlon_earth(t::Number)
  	mod2pi(@evalpoly(t, 1.753470314, 628.3075849991))
end


"""
	fa_mlon_mars(t::Number) 

Return the mean heliocentric longitude of Mars in radians, at time `t` expressed in `TDB` 
Julian centuries since [`J2000`](@ref). 

### References 
- Luzum, B. and Petit G. (2012). The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/fama03.c) software library
"""
function fa_mlon_mars(t::Number)
	mod2pi(@evalpoly(t, 6.203480913, 334.0612426700))
end


"""
	fa_mlon_jupiter(t::Number) 

Return the mean heliocentric longitude of Jupiter in radians, at time `t` expressed in `TDB`
Julian centuries since [`J2000`](@ref). 

### References 
- Luzum, B. and Petit G. (2012). The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/faju03.c) software library
"""
function fa_mlon_jupiter(t::Number)
  	mod2pi(@evalpoly(t, 0.599546497, 52.9690962641))
end


"""
	fa_mlon_saturn(t::Number) 

Return the mean heliocentric longitude of Saturn in radians, at time `t` expressed in `TDB` 
Julian centuries since [`J2000`](@ref). 

### References 
- Luzum, B. and Petit G. (2012). The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/fasa03.c) software library
"""
function fa_mlon_saturn(t::Number)
  	mod2pi(@evalpoly(t, 0.874016757, 21.3299104960))
end


"""
	fa_mlon_uranus(t::Number) 

Return the mean heliocentric longitude of Uranus in radians, at time `t` expressed in `TDB` 
Julian centuries since [`J2000`](@ref). 

### References 
- Luzum, B. and Petit G. (2012). The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/faur03.c) software library
"""
function fa_mlon_uranus(t::Number)
  	mod2pi(@evalpoly(t, 5.481293872, 7.4781598567))
end


"""
	fa_mlon_neptune(t::Number) 

Return the mean heliocentric longitude of Neptune in radians, at time `t` expressed in `TDB` 
Julian centuries since [`J2000`](@ref). 

### References 
- Luzum, B. and Petit G. (2012). The IERS Conventions (2010), 
  [IERS Technical Note No. 36](https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html) 
- [ERFA](https://github.com/liberfa/erfa/blob/master/src/fane03.c) software library
"""
function fa_mlon_neptune(t::Number)
  	mod2pi(@evalpoly(t, 5.311886287, 3.8133035638))
end


