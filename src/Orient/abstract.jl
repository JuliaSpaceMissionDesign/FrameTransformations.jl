import Basic.Bodies: body_naifid, CelestialBody, NAIFId

"""
    orient_declination_rate(body::NAIFId, T::Float64)

Return the right ascension of the body-fixed frame of `body` w.r.t. the ICRF 
at epoch `ep` in rad. Here `T` is expressed in centuries since J2000 in the 
TDB time scale.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)
"""
function orient_right_ascension end

"""
    orient_declination_rate(body::NAIFId, T::Float64)

Return the right ascension rate of change of the body-fixed frame of `body` w.r.t. the
ICRF at epoch `ep` in rad/s. Here `T` is expressed in centuries since J2000 in the 
TDB time scale.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)
"""
function orient_right_ascension_rate end

"""
    orient_declination_rate(body::NAIFId, T::Float64)

Return the declination of the body-fixed frame of `body` w.r.t. the ICRF 
at epoch `ep` in rad. Here `T` is expressed in centuries since J2000 in the 
TDB time scale.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)
"""
function orient_declination end

"""
    orient_declination_rate(body::NAIFId, T::Float64)

Return the declination rate of change of the body-fixed frame of `body` w.r.t. the
ICRF at epoch `ep` in rad/s. Here `T` is expressed in centuries since J2000 in the 
TDB time scale.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)
"""
function orient_declination_rate end

"""
    orient_rotation_angle(body::NAIFId, T::Float64)

Return the declination of the body-fixed frame of `body` w.r.t. the ICRF 
at epoch `ep` in rad. Here `T` is expressed in centuries since J2000 in the 
TDB time scale.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)
"""
function orient_rotation_angle end

"""
    orient_rotation_rate(body::NAIFId, T::Float64)

Return the declination rate of change of the body-fixed frame of `body` w.r.t. the
ICRF at epoch `ep` in rad/s. Here `T` is expressed in centuries since J2000 in the 
TDB time scale.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)
"""
function orient_rotation_rate end

for fun in (:orient_declination, :orient_declination_rate, 
    :orient_right_ascension, :orient_right_ascension_rate, 
    :orient_rotation_angle, :orient_rotation_rate)
    @eval begin
        $fun(id::Integer, T::N) where {N<:AbstractFloat} = $fun(NAIFId(id), T)
        $fun(naifid::NAIFId, T::N) where {N<:AbstractFloat} = $fun(naifid.valId, T)
        $fun(body::B, T::N) where {B<:CelestialBody, N<:AbstractFloat} = $fun(body_naifid(body), T)
        export $fun
    end
end