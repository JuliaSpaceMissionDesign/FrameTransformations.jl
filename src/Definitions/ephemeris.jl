
@interface function add_point_ephemeris!(
    ::FrameSystem{O, N}, ::AbstractEphemerisProvider, ::Symbol, ::Int
) where {O, N} end

@interface function add_point_ephemeris!(
    ::FrameSystem{O, N}, ::AbstractEphemerisProvider, ::Dict{Int, Symbol}
) where {O, N} end

