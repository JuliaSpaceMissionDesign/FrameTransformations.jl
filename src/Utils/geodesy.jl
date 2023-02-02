"""
    geoc2pos(h::Number, λ::Number, ϕ::Number, R::Number)

Transform geocentric coordinates in a cartesian position vector, given the longitude `λ`, 
the geocentric latitude `ϕ`, the geocentric altitude over the reference sphere `h` and its 
radius `R`.
"""
@fastmath function geoc2pos(h::Number, λ::Number, ϕ::Number, R::Number)
    sϕ, cϕ = sincos(ϕ)
    sλ, cλ = sincos(λ)
    r = R + h
    return SA[r*cϕ*cλ, r*cϕ*sλ, r*sϕ]
end


"""
    pos2geoc(pos::AbstractVector, R::Number)

Transform a cartesian 3-elements position vector `pos` into longitude, geocentric latitude
and altitude over the reference sphere with radius `R`, respectively.
"""
@fastmath function pos2geoc(pos::AbstractVector, R::Number)
    @inbounds begin 
        r = sqrt(pos[1]^2 + pos[2]^2 + pos[3]^2)

        ϕ = asin(pos[3]/r)
        λ = atan(pos[2], pos[1])
        h = r - R
    end
    return λ, ϕ, h
end


"""
    pos2geod(pos::AbstratVector, R::Number, f::Number, toll::Number=1e-12)

Transform a cartesian 3-elements position vector `pos` into longitude, geodetic latitude 
and altitude over the reference ellipsoid with radius `R` and flattening `f`. 
 
### References 
- Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications, Microcosm Press, 
    Hawthorn, CA, USA.
"""
@fastmath function pos2geod(pos::AbstractVector, R::Number, f::Number, toll::Number=1e-12)
    
    @inbounds x, y, z = pos[1], pos[2], pos[3]
    sz = sign(z)

    rδs = x^2 + y^2

    # Get eccentricity from flattening  
    e² = f*(2-f)
    sϕ² = z^2/(rδs + z^2) 
    
    zₙ, cₙ = z, 1.0
    err = toll + 1
    while err > toll
        c = cₙ

        z² = zₙ^2 
        sϕ² = z²/(rδs + z²)

        cₙ = R*e²*sqrt(sϕ²/(1-e²*sϕ²))  
        zₙ = z + cₙ*sz 

        err = abs(cₙ - c)
    end

    λ = atan(y, x)
    ϕ = atan(zₙ, sqrt(rδs))
    hₑ = sqrt(rδs + zₙ^2) - R/sqrt(1-e²*sϕ²)
    
    return λ, ϕ, hₑ

end


"""
    geod2pos(h::Number, λ::Number, ϕ::Number, R::Number, f::Number)

Transform longitude `λ`, geodetic latitude `ϕ` and altitude over the reference ellipsoid to 
a cartesian position vector, given the reference radius `R` and the flattening `f`.

### References 
- Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications, Microcosm Press, 
    Hawthorn, CA, USA.
""" 
@fastmath function geod2pos(h::Number, λ::Number, ϕ::Number, R::Number, f::Number)

    # Get eccentricity from flattening 
    e² = (2-f)*f

    sϕ, cϕ = sincos(ϕ)
    sλ, cλ = sincos(λ)

    d = R/sqrt(1-e²*sϕ^2)
    c = (d+h)*cϕ
    s = (1-e²)*d

    SA[c*cλ, c*sλ, (s+h)*sϕ]
end
