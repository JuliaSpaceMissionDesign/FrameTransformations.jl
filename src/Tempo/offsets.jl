using ..Utils

function Epoch(ep::Epoch{S1}, ::S2) where {S1<:TimeScale, S2<:TimeScale}
    second, fraction, error = apply_offset(ep.second, ep.fraction, ep.error, S1(), S2())
    Epoch{S2}(second, fraction, error)
end

@inline function Utils.apply_offset(second::Int64,
    fraction::T,
    error::T,
    from::S1,
    to::S2)::Tuple{Int64, T, T} where {T, S1<:TimeScale, S2<:TimeScale}
    path = find_path(from, to)
    isempty(path) && throw(NoPathError(string(from), string(to)))
    length(path) == 2 && return _apply_offset(second, fraction, error, from, to)
    return _apply_offset((second, fraction, error), path...)
end

@inline function _apply_offset(second::Int64,
    fraction::T,
    error::T,
    from::S1,
    to::S2)::Tuple{Int64, T, T} where {T, S1<:TimeScale, S2<:TimeScale}
    return Utils.apply_offset(second, fraction, error, getoffset(from, to, second, fraction))
end

@generated function _apply_offset(sf, path...)
    expr = :(sf)
    for i in 1:length(path) - 1
        s1 = path[i]
        s2 = path[i+1]
        expr = :(_apply_offset($expr..., $s1(), $s2()))
    end
    return quote
        Base.@_inline_meta
        $expr
    end
end

######
# TT #
######

const OFFSET_TAI_TT = 32.184  # seconds 

"""
    getoffset(TAI, TT, args...)
Return the fixed offset between [`TAI`](@ref) and [`TT`](@ref) in seconds.
"""
getoffset(::InternationalAtomicTime, ::TerrestrialTime, args...) = OFFSET_TAI_TT

"""
    getoffset(TT, TAI, args...)
Return the fixed offset between [`TT`](@ref) and [`TAI`](@ref) in seconds.
"""
getoffset(::TerrestrialTime, ::InternationalAtomicTime, args...) = -OFFSET_TAI_TT

#######
# TCG #
#######

const JD77_SEC = -7.25803167816e8
const LG_RATE = 6.969290134e-10

"""
    getoffset(TCG, TT, second, fraction)
Return the linear offset between [`TCG`](@ref) and [`TT`](@ref) for the
current epoch (`second` after J2000 and `fraction`) in seconds.
"""
function getoffset(::GeocentricCoordinateTime, ::TerrestrialTime, second, fraction)
    dt = second - JD77_SEC + fraction
    return -LG_RATE * dt
end

"""
    getoffset(TT, TCG, second, fraction)
Return the linear offset between [`TT`](@ref) and [`TCG`](@ref) for the
current epoch (`second` after J2000 and `fraction`) in seconds.
"""
function getoffset(::TerrestrialTime, ::GeocentricCoordinateTime, second, fraction)
    rate = LG_RATE / (1.0 - LG_RATE)
    dt = second - JD77_SEC + fraction
    return rate * dt
end

#######
# TCB #
#######

const LB_RATE = 1.550519768e-8

"""
    getoffset(TCB, TDB, second, fraction)
Return the linear offset between [`TCB`](@ref) and [`TDB`](@ref) for the
current epoch (`second` after J2000 and `fraction`) in seconds.
"""
function getoffset(::BarycentricCoordinateTime, ::BarycentricDynamicalTime, second, fraction)
    dt = second - JD77_SEC + fraction
    return -LB_RATE * dt
end

"""
    getoffset(TDB, TCB, second, fraction)
Return the linear offset between [`TDB`](@ref) and [`TCB`](@ref) for the
current epoch (`second` after J2000 and `fraction`) in seconds.
"""
function getoffset(::BarycentricDynamicalTime, ::BarycentricCoordinateTime, second, fraction)
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
    getoffset(TT, TDB, second, fraction[, eop])
Return the offset between [`TT`](@ref) and [`TDB`](@ref) for the
current epoch (`second` after J2000 and `fraction`) in seconds.
This routine is accurate to ~40 microseconds over the interval 1900-2100.
!!! note
    An accurate transformation between TDB and TT depends on the
    trajectory of the observer. For two observers fixed on Earth's surface
    the quantity TDB-TT can differ by as much as ~4 microseconds. 
# Example
### References ###
- [https://www.cv.nrao.edu/~rfisher/Ephemerides/times.html#TDB](https://www.cv.nrao.edu/~rfisher/Ephemerides/times.html#TDB)
- [Issue #26](https://github.com/JuliaAstro/AstroTime.jl/issues/26)
"""
@inline function getoffset(::TerrestrialTime, ::BarycentricDynamicalTime,
                           second, fraction)
    tt = fraction + second
    g = m₀ + m₁ * tt
    return k * sin(g + eb * sin(g))
end

"""
    getoffset(TDB, TT, second, fraction[, eop])
Return the offset between [`TDB`](@ref) and [`TT`](@ref) for the
current epoch (`second` after J2000 and `fraction`) in seconds.
This routine is accurate to ~40 microseconds over the interval 1900-2100.
!!! note
    An accurate transformation between TDB and TT depends on the
    trajectory of the observer. For two observers fixed on Earth's surface
    the quantity TDB-TT can differ by as much as ~4 microseconds. 
# Example
### References ###
- [https://www.cv.nrao.edu/~rfisher/Ephemerides/times.html#TDB](https://www.cv.nrao.edu/~rfisher/Ephemerides/times.html#TDB)
- [Issue #26](https://github.com/JuliaAstro/AstroTime.jl/issues/26)
"""
@inline function getoffset(::BarycentricDynamicalTime, ::TerrestrialTime,
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
