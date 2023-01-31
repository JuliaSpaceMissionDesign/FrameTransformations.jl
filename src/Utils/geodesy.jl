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
"""
@fastmath function pos2geod(pos::AbstractVector, R::Number, f::Number, toll::Number=1e-12)
    
    # Get eccentricity from flattening 
    e = sqrt(1-(1-f)^2)
    e² = e*e 

    @inbounds ri, rj, rk = pos[1], pos[2], pos[3]

    rδs = ri^2 + rj^2
    r = rδs + rk^2 
    rδ = sqrt(rδs)

    sϕ = rk/r 
    ϕ = asin(sϕ)

    d = sqrt(1-e²*sϕ^2)
    c = R/d

    err, iter = 1.0, 1
    while iter < 5 && err > toll
        ϕn = atan(rk + c*e²*sϕ, rδ)
        err, ϕ = abs(ϕn - ϕ), ϕn

        sϕ = sin(ϕ)
        d = sqrt(1-e²*sϕ^2)
        c = R/d
        
        iter += 1
    end

    if π/2 - abs(ϕ) < π/180
        s = R*(1-e²)/d
        hₑ = rk/sϕ - s
    else 
        cϕ = cos(ϕ)
        hₑ = rδ/cϕ - c 
    end
    
    λ = atan(rj, ri)
    return λ, ϕ, hₑ

end


"""
    geod2pos(h::Number, λ::Number, ϕ::Number, R::Number, f::Number)

Transform longitude `λ`, geodetic latitude `ϕ` and altitude over the reference ellipsoid to 
a cartesian position vector, given the reference radius `R` and the flattening `f`.
""" 
@fastmath function geod2pos(h::Number, λ::Number, ϕ::Number, R::Number, f::Number)

    # Get eccentricity from flattening 
    e = sqrt(1-(1-f)^2)
    e² = e*e 

    sϕ, cϕ = sincos(ϕ)
    sλ, cλ = sincos(λ)

    d = sqrt(1-e²*sϕ^2)
    c, s = R/d, R*(1-e²)/d 

    SA[(c+h)*cϕ*cλ, (c+h)*cϕ*sλ, (s+h)*sϕ]
end
