"""
    AbstractEphemerisProvider

Abstract type to represent ephemeris types.
"""
abstract type AbstractEphemerisProvider end

"""
    ephem_load(::Type{E}, file::String) where {E<:AbstractEphemerisProvider}

Load a generic ephemeris file.
"""
function ephem_load(::Type{E}, file::String) where {E<:AbstractEphemerisProvider}
    return E(file)
end


for fun in (:ephem_position_records, 
            :ephem_orient_records, 
            :ephem_available_points, 
            :ephem_available_axes, 
            :ephem_timespan, 
            :ephem_timescale, 
            :ephem_compute_order!, 
            :ephem_orient_order!
            )

    @eval begin 
        function ($fun)(E::AbstractEphemerisProvider)
            throw(NotImplementedError(
                String(Symbol(@__MODULE__)),
                "$($fun) shall be implemented for $(typeof(E))"
            ))
        end
    end

end 
