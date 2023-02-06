"""
    geoc2pos(r::Number, λ::Number, ϕ::Number)
    geoc2pos(geoc::AbstractArray)

Transform geocentric coordinates in a cartesian position vector, given the longitude `λ`, 
the geocentric latitude `ϕ` and the radius `r`.
"""
@fastmath function geoc2pos(r::Number, λ::Number, ϕ::Number)
    sϕ, cϕ = sincos(ϕ)
    sλ, cλ = sincos(λ)
    return SA[r*cϕ*cλ, r*cϕ*sλ, r*sϕ]
end

@inline geoc2pos(geoc::AbstractVector) = geoc2pos(geoc[1], geoc[2], geoc[3])

"""
    pos2geoc(pos::AbstractVector)

Transform a cartesian 3-elements position vector `pos` into radius, longitude and geocentric 
latitude, respectively.
"""
@fastmath function pos2geoc(pos::AbstractVector)

    @inbounds x, y, z = pos[1], pos[2], pos[3]
    r = sqrt(x^2 + y^2 + z^2)

    λ = atan(y, x)
    ϕ = asin(z/r)
    
    return SA[r, λ, ϕ]
end


"""
    geoc2pv(geoc::AbstractVector)

Transform a spherical geocentric 6-elements state vector (radius, longitude, geocentric 
latitude and their derivatives) into a cartesian 6-elements vector (position and velocity).
"""
function geoc2pv(geoc::AbstractVector)

    @fastmath @inbounds begin 
        r = geoc[1]   
        δr, δλ, δϕ = geoc[4], geoc[5], geoc[6]

        sλ, cλ = sincos(geoc[2])
        sϕ, cϕ = sincos(geoc[3])
    end

    δϕr = δϕ*r
    δϕrsϕ = δϕr*sϕ

    δrcϕ = δr*cϕ
    rcϕδλ = r*cϕ*δλ

    SA[r*cϕ*cλ, r*cϕ*sλ, r*sϕ, 
        (δrcϕ - δϕrsϕ)*cλ - rcϕδλ*sλ, 
        (δrcϕ - δϕrsϕ)*sλ + rcϕδλ*cλ, 
        δr*sϕ + δϕr*cϕ]

end

"""
    pv2geoc(pv::AbstractVector)

Transform a cartesian 6-elements state vector (position and velocity) into radius, longitude, 
geocentric latitude and their derivatives, respectively.
"""
@fastmath function pv2geoc(pv::AbstractVector)
    @inbounds x, y, z = pv[1], pv[2], pv[3]
    @inbounds dx, dy, dz = pv[4], pv[5], pv[6]

    rxy2 = x^2 + y^2
    rxy = sqrt(rxy2)

    r2 = rxy2 + z^2
    r = sqrt(r2)

    # radius
    δr = (x*dx + y*dy + z*dz)/r
    δrxy = (x*dx + y*dy)/rxy

    # longitude
    λ = atan(y, x)
    δλ = (dy*x - dx*y)/rxy2

    # latitude
    ϕ = atan(z, rxy)
    δϕ = (dz*rxy - z*δrxy)/r2

    return SA[r, λ, ϕ, δr, δλ, δϕ]
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
