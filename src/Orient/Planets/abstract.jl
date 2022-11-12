import Basic.Bodies: CelestialBody

for fun in (
    :orient_right_ascension, :orient_right_ascension_rate, 
    :orient_declination, :orient_declination_rate,
    :orient_rotation_angle, :orient_rotation_rate)
    @eval begin
        function $fun(b::CelestialBody, TDB::Number) 
            throw(
                NotImplementedError(
                    String(Symbol(@__MODULE__)), 
                    "`$($(fun))` shall be implemented for $(b)"
                )
            )
        end
    end
end

#
#   Documentation
#

"""
    orient_right_ascension(::CelestialBody, TDB::Number)

Return the right ascension of the body-fixed frame of `body` w.r.t. the ICRF 
at epoch `ep` in rad. Here `TDB` is expressed in centuries since J2000 in the 
TDB time scale.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)

!!! warning 
    This method is abstract! A concrete implementation for the desired body 
    shall be defined.
"""
orient_right_ascension

"""
    orient_right_ascension_rate(::CelestialBody, TDB::Number)

Return the right ascension rate of change of the body-fixed frame of `body` w.r.t. the
ICRF at epoch `ep` in rad/s. Here `TDB` is expressed in centuries since J2000 in the 
TDB time scale.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)

!!! warning 
    This method is abstract! A concrete implementation for the desired body 
    shall be defined.
"""
orient_right_ascension_rate

"""
    orient_declination(::CelestialBody, TDB::Number)

Return the declination of the body-fixed frame of `body` w.r.t. the ICRF 
at epoch `ep` in rad. Here `TDB` is expressed in centuries since J2000 in the 
TDB time scale.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)

!!! warning 
    This method is abstract! A concrete implementation for the desired body 
    shall be defined.
"""
orient_declination

"""
    orient_declination_rate(::CelestialBody, TDB::Number)

Return the declination rate of change of the body-fixed frame of `body` w.r.t. the
ICRF at epoch `ep` in rad/s. Here `TDB` is expressed in centuries since J2000 in the 
TDB time scale.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)

!!! warning 
    This method is abstract! A concrete implementation for the desired body 
    shall be defined.
"""
orient_declination_rate

"""
    orient_rotation_angle(::CelestialBody, TDB::Number)

Return the declination of the body-fixed frame of `body` w.r.t. the ICRF 
at epoch `ep` in rad. Here `TDB` is expressed in centuries since J2000 in the 
TDB time scale.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)

!!! warning 
    This method is abstract! A concrete implementation for the desired body 
    shall be defined.
"""
orient_rotation_angle

"""
    orient_rotation_rate(::CelestialBody, TDB::Number)

Return the declination rate of change of the body-fixed frame of `body` w.r.t. the
ICRF at epoch `ep` in rad/s. Here `TDB` is expressed in centuries since J2000 in the 
TDB time scale.

# References
- Archinal, Brent Allen, et al. "Report of the IAU Working Group on Cartographic Coordinates 
   and Rotational Elements: 2015." *Celestial Mechanics and Dynamical Astronomy* 
   volume 130, Article number: 22 (2018)

!!! warning 
    This method is abstract! A concrete implementation for the desired body 
    shall be defined.
"""
orient_rotation_rate