import Basic: NotImplementedError

for fun in (
    :body_naifid, :body_parent, :body_system_equivalent, :body_gm,
    :body_mean_radius, :body_polar_radius, :body_equatorial_radius,
    :body_flattening)
    @eval begin
        function $fun(b::CelestialBody) 
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
# Documentation
#

"""
    body_naifid(body::CelestialBody)::Integer

Return the NAIF ID code for `body`. 

!!! warning 
    This method is abstract! A concrete implementation for the desired body 
    shall be defined.
"""
body_naifid

"""
    body_parent(body::CelestialBody)::Integer

Get parent of a given body.
Returns an integer, which is the `astroid` of the parent body.

!!! warning 
    This method is abstract! A concrete implementation for the desired body 
    shall be defined.
"""
body_parent

"""
    body_system_equivalent(body::CelestialBody)::Float64

Return the body system equivalent body or barycenter.

!!! warning 
    This method is abstract! A concrete implementation for the desired body 
    shall be defined.
"""
body_system_equivalent

"""
    body_gm(body::CelestialBody)::Float64

Return the gravitational parameter ``\\mu = GM`` of `body` in km^3/s^2.

!!! warning 
    This method is abstract! A concrete implementation for the desired body 
    shall be defined.
"""
body_gm

"""
    body_mean_radius(body::CelestialBody)::Float64

Return the mean radius of `body` in km.

!!! warning 
    This method is abstract! A concrete implementation for the desired body 
    shall be defined.
"""
body_mean_radius

"""
    body_polar_radius(body::CelestialBody)::Float64

Return the polar radius of `body` in km.

!!! warning 
    This method is abstract! A concrete implementation for the desired body 
    shall be defined.
"""
body_polar_radius

"""
    body_equatorial_radius(body::CelestialBody)::Float64

Return the polar radius of `body` in km.

!!! warning 
    This method is abstract! A concrete implementation for the desired body 
    shall be defined.
"""
body_equatorial_radius

"""
    body_flattening(body::CelestialBody)::Float64

Return the flattening of `body`.

!!! warning 
    This method is abstract! A concrete implementation for the desired body 
    shall be defined.
"""
body_flattening