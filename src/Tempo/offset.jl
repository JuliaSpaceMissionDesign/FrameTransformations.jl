######
# TT #
######

const OFFSET_TAI_TT = 32.184  # seconds 

"""
    offset_tai2tt(seconds)

Return the fixed offset between [`TAI`](@ref) and [`TT`](@ref) in seconds.
"""
@inline offset_tai2tt(seconds) = OFFSET_TAI_TT

"""
    offset_tt2tai(seconds)

Return the fixed offset between [`TT`](@ref) and [`TAI`](@ref) in seconds.
"""
@inline offset_tt2tai(seconds) = -OFFSET_TAI_TT

#######
# TCG #
#######

const JD77_SEC = -7.25803167816e8
const LG_RATE = 6.969290134e-10

"""
    offset_tcg2tt(seconds)

Return the linear offset between [`TCG`](@ref) and [`TT`](@ref) in seconds.
"""
@inline function offset_tcg2tt(seconds)
    δt = seconds - JD77_SEC
    return -LG_RATE * δt
end

"""
    offset_tt2tcg(seconds)

Return the linear offset between [`TT`](@ref) and [`TCG`](@ref) in seconds.
"""
@inline function offset_tt2tcg(seconds)
    rate = LG_RATE / (1.0 - LG_RATE)
    δt = seconds - JD77_SEC 
    return rate * δt
end

#######
# TCB #
#######

const LB_RATE = 1.550519768e-8

"""
    offset_tcb2tdb(seconds)
    
Return the linear offset between [`TCB`](@ref) and [`TDB`](@ref) in seconds.
"""
@inline function offset_tcb2tdb(seconds)
    δt = seconds - JD77_SEC
    return -LB_RATE * δt
end

"""
    offset_tdb2tcb(seconds)
    
Return the linear offset between [`TDB`](@ref) and [`TCB`](@ref) in seconds.
"""
@inline function offset_tdb2tcb(seconds)
    rate = LB_RATE / (1.0 - LB_RATE)
    δt = seconds - JD77_SEC
    return rate * δt
end

#######
# TDB #
#######

const k = 1.657e-3
const eb = 1.671e-2
const m₀ = 6.239996
const m₁ = 1.99096871e-7

"""
    offset_tt2tdb(seconds)

Return the offset between [`TT`](@ref) and [`TDB`](@ref) in seconds.
This routine is accurate to ~40 microseconds over the interval 1900-2100.

!!! note
    An accurate transformation between TDB and TT depends on the trajectory of the observer. 
    For two observers fixed on Earth's surface the quantity TDB-TT can differ by as much 
    as ~4 microseconds.

### References 
- https://www.cv.nrao.edu/~rfisher/Ephemerides/times.html#TDB
- [Issue #26](https://github.com/JuliaAstro/AstroTime.jl/issues/26)
"""
@inline function offset_tt2tdb(seconds)
    g = m₀ + m₁ * seconds
    return k * sin(g + eb * sin(g))
end

"""
    offset_tdb2tt(seconds)

Return the offset between [`TDB`](@ref) and [`TT`](@ref) in seconds.
This routine is accurate to ~40 microseconds over the interval 1900-2100.

!!! note
    An accurate transformation between TDB and TT depends on the trajectory of the observer. 
    For two observers fixed on Earth's surface the quantity TDB-TT can differ by as much 
    as ~4 microseconds.

### References 
- https://www.cv.nrao.edu/~rfisher/Ephemerides/times.html#TDB
- [Issue #26](https://github.com/JuliaAstro/AstroTime.jl/issues/26)
"""
@inline function offset_tdb2tt(seconds)
    tt = seconds
    offset = 0.0
    for _ in 1:3
        g = m₀ + m₁ * tt
        offset = -k * sin(g + eb * sin(g))
        tt = seconds + offset
    end
    return offset
end

#######
# UTC #
#######

"""
    offset_tai2utc(seconds)

Return the offset between [`TAI`](@ref) and [`UTC`](@ref) in seconds.
"""
@inline function offset_tai2utc(seconds)
    tai = seconds/86400.0
    _, utc = tai2utc(DJ2000, tai)
    return (utc - tai) * 86400
end

"""
    offset_utc2tai(seconds)

Return the offset between [`UTC`](@ref) and [`TAI`](@ref) in seconds.
"""
@inline function offset_utc2tai(seconds)
    utc = seconds/86400.0 
    _, tai = utc2tai(DJ2000, utc)
    return (tai - utc) * 86400.0
end

#######
# UT1 #
#######

# using Basic.Orient: IERS_EOP

@inline function offset_utc2ut1(seconds)
    utc = seconds/86400.0
    return IERS_EOP.UT1_UTC[utc]
end