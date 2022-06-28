"""
    AbstractEphemeris

Abstract type to represent ephemeris types.
"""
abstract type AbstractEphemeris end


"""
    vector3!(r::Vector{NV}, eph::AbstractEphemeris, from::Int, to::Int, jd2000::N, [time]::N) where {N, NV}

Compute position vector from the center `from` to the target `to`. 

### Input/s

- `r` -- Position vector (**output**)
- `eph` -- Ephemeris object type. This must be a concrete implementation
           of `AbstractEphemeris`
- `from` -- observer
- `to` -- target
- `jd2000` -- integer part of the date (in TDB scale)
- `time` -- fractional part of the date (in TDB scale), optional
"""
function vector3! end 

"""
    vector6!(rv::Vector{NV}, eph::AbstractEphemeris, from::Int, to::Int, jd2000::N, [time]::N) where {N, NV}

Compute position and velocity vector from the center `from` to the target `to`. 

### Input/s

- `rv` -- Position and velocity vector (**output**)
- `eph` -- Ephemeris object type. This must be a concrete implementation
           of `AbstractEphemeris`
- `from` -- observer
- `to` -- target
- `jd2000` -- integer part of the date (in TDB scale)
- `time` -- fractional part of the date (in TDB scale), optional
"""
function vector6! end 

"""
    vector9!(rva::Vector{NV}, eph::AbstractEphemeris, from::Int, to::Int, jd2000::N, [time]::N) where {N, NV}

Compute position, velocity and acceleration vector from the center 
`from` to the target `to`. 

### Input/s

- `rva` -- Position, velocity and acceleration vector (**output**)
- `eph` -- Ephemeris object type. This must be a concrete implementation
           of `AbstractEphemeris`
- `from` -- observer
- `to` -- target
- `jd2000` -- integer part of the date (in TDB scale)
- `time` -- fractional part of the date (in TDB scale), optional
"""
function vector9! end 

"""
    vector12!(rva::Vector{NV}, eph::AbstractEphemeris, from::Int, to::Int, jd2000::N, [time]::N) where {N, NV}

Compute position, velocity, acceleration and jerk vector from the center 
`from` to the target `to`. 

### Input/s

- `rvaj` -- Position, velocity, acceleration and jerk vector (**output**)
- `eph` -- Ephemeris object type. This must be a concrete implementation
           of `AbstractEphemeris`
- `from` -- observer
- `to` -- target
- `jd2000` -- integer part of the date (in TDB scale)
- `time` -- fractional part of the date (in TDB scale), optional
"""
function vector12! end 

function vector3 end 
function vector6 end 
function vector9 end 
function vector12 end 
