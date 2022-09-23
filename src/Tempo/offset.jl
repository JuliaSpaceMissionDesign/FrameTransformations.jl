@inline function apply_offset(s::Int64, f, s2f)
    # accurate at 1e-12 seconds
    s2 = floor(Int64, s2f)
    f2 = s2f - s2
    sec = s + s2 
    fi = floor(Int64, f)
    f2i = floor(Int64, f2)

    ff = f - fi 
    f2f = f2 - f2i 
    fri = floor(Int64, ff + f2f)

    return sec + fi + f2i + fri, ff + f2f - fri
end

using StaticArrays


"""
    apply_offsets(second::N, fraction::T, 
    from::S1, to::S2)::Tuple{N, T} where {N<:Integer,
    S1, S2, T<:AbstractFloat}

Apply offset to transform `seconds` and `fraction` since J2000 from time scale 
`S1` to time scale `S1`.
"""
@inline function apply_offsets(second::N, fraction::T,
    from::S1, to::S2)::Tuple{N, T} where {N<:Integer,
    S1, S2, T<:AbstractFloat}

    path = find_path(from, to)
    if isempty(path) 
        throw(error("[Tempo] no path between scales $from and $to"))
    end
    if length(path) == 2
        return _apply_offset(second, fraction, from, to)
    end
    @inbounds for i in 1:length(path)-1
        s1 = path[i]
        s2 = path[i+1]
        second, fraction = _apply_offset(second, fraction, s1, s2)
    end
    return second, fraction
end

@inline function _apply_offset(second::N, fraction::T, 
    from::S1, to::S2)::Tuple{N, T} where {T, N<:Integer, S1, S2}
    return apply_offset(second, fraction, offset(from, to, second, fraction))
end

######
# TT #
######

const OFFSET_TAI_TT = 32.184  # seconds 

"""
    offset(TAI, TT, args...)

Return the fixed offset between [`TAI`](@ref) and [`TT`](@ref) in seconds.
"""
offset(::InternationalAtomicTime, ::TerrestrialTime, second, fraction) = OFFSET_TAI_TT

"""
    offset(TT, TAI, args...)

Return the fixed offset between [`TT`](@ref) and [`TAI`](@ref) in seconds.
"""
offset(::TerrestrialTime, ::InternationalAtomicTime, second, fraction) = -OFFSET_TAI_TT

#######
# TCG #
#######

const JD77_SEC = -7.25803167816e8
const LG_RATE = 6.969290134e-10

"""
    offset(TCG, TT, second, fraction)

Return the linear offset between [`TCG`](@ref) and [`TT`](@ref) for the
current epoch (`second` after J2000 and `fraction`) in seconds.
"""
function offset(::GeocentricCoordinateTime, ::TerrestrialTime, second, fraction)
    dt = second - JD77_SEC + fraction
    return -LG_RATE * dt
end

"""
    offset(TT, TCG, second, fraction)

Return the linear offset between [`TT`](@ref) and [`TCG`](@ref) for the
current epoch (`second` after J2000 and `fraction`) in seconds.
"""
function offset(::TerrestrialTime, ::GeocentricCoordinateTime, second, fraction)
    rate = LG_RATE / (1.0 - LG_RATE)
    dt = second - JD77_SEC + fraction
    return rate * dt
end

#######
# TCB #
#######

const LB_RATE = 1.550519768e-8

"""
    offset(TCB, TDB, second, fraction)
    
Return the linear offset between [`TCB`](@ref) and [`TDB`](@ref) for the
current epoch (`second` after J2000 and `fraction`) in seconds.
"""
function offset(::BarycentricCoordinateTime, ::BarycentricDynamicalTime, second, fraction)
    dt = second - JD77_SEC + fraction
    return -LB_RATE * dt
end

"""
    offset(TDB, TCB, second, fraction)

Return the linear offset between [`TDB`](@ref) and [`TCB`](@ref) for the
current epoch (`second` after J2000 and `fraction`) in seconds.
"""
function offset(::BarycentricDynamicalTime, ::BarycentricCoordinateTime, second, fraction)
    rate = LB_RATE / (1.0 - LB_RATE)
    dt = second - JD77_SEC + fraction
    return rate * dt
end

#######
# TDB #
#######

const k = 1.657e-3
const eb = 1.671e-2
const m₀ = 6.239996
const m₁ = 1.99096871e-7

"""
    offset(TT, TDB, second, fraction[, eop])

Return the offset between [`TT`](@ref) and [`TDB`](@ref) for the
current epoch (`second` after J2000 and `fraction`) in seconds.
This routine is accurate to ~40 microseconds over the interval 1900-2100.

!!! note
    An accurate transformation between TDB and TT depends on the
    trajectory of the observer. For two observers fixed on Earth's surface
    the quantity TDB-TT can differ by as much as ~4 microseconds.


### References 
- [https://www.cv.nrao.edu/~rfisher/Ephemerides/times.html#TDB](https://www.cv.nrao.edu/~rfisher/Ephemerides/times.html#TDB)
- [Issue #26](https://github.com/JuliaAstro/AstroTime.jl/issues/26)
"""
@inline function offset(::TerrestrialTime, ::BarycentricDynamicalTime,
                           second, fraction)
    tt = fraction + second
    g = m₀ + m₁ * tt
    return k * sin(g + eb * sin(g))
end

"""
    offset(TDB, TT, second, fraction[, eop])

Return the offset between [`TDB`](@ref) and [`TT`](@ref) for the
current epoch (`second` after J2000 and `fraction`) in seconds.
This routine is accurate to ~40 microseconds over the interval 1900-2100.

!!! note
    An accurate transformation between TDB and TT depends on the
    trajectory of the observer. For two observers fixed on Earth's surface
    the quantity TDB-TT can differ by as much as ~4 microseconds. 

### References
- [https://www.cv.nrao.edu/~rfisher/Ephemerides/times.html#TDB](https://www.cv.nrao.edu/~rfisher/Ephemerides/times.html#TDB)
- [Issue #26](https://github.com/JuliaAstro/AstroTime.jl/issues/26)
"""
@inline function offset(::BarycentricDynamicalTime, ::TerrestrialTime,
                           second, fraction)
    tdb = fraction + second
    tt = tdb
    offset = 0.0
    for _ in 1:3
        g = m₀ + m₁ * tt
        offset = -k * sin(g + eb * sin(g))
        tt = tdb + offset
    end
    return offset
end

#######
# UTC #
#######

@inline function offset(::InternationalAtomicTime, ::CoordinatedUniversalTime, 
                           seconds, fraction)
    tai = (fraction + seconds)/86400.0
    _, utc = tai2utc(DJ2000, tai)
    return (utc - tai)*86400.0
end
   
@inline function offset(::CoordinatedUniversalTime, ::InternationalAtomicTime,
                           seconds, fraction)
    utc = (fraction + seconds)/86400.0
    _, tai = utc2tai(DJ2000, utc)
    return (tai - utc)*86400.0
end