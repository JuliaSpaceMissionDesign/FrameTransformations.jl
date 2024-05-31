module CalcephEphemerisExt

import FrameTransformations: add_point_ephemeris!
using FrameTransformations: FrameSystem, add_point!, check_point_ephemeris, 
                            FramePointFunctions, POINT_CLASSID_DYNAMIC
using JSMDInterfaces.Ephemeris: ephem_compute!    
using CalcephEphemeris: CalcephProvider

function add_point_ephemeris!(
    frames::FrameSystem{O, N}, eph::CalcephProvider, name::Symbol, id::Int
) where {O, N}

    parentid, axesid = check_point_ephemeris(frames, eph, id)
    
    # Note: this cache would not work with ForwardDiff AD types but this is intentional as 
    # Calceph would be not compatible
    nth = Threads.nthreads()
    cache = [zeros(N, 3O) for _ in 1:nth]

    funs = FramePointFunctions{O, N}(
        t -> SVectorNT{3O, N}( _wrap_vector3(cache, eph, parentid, id, t) ),
        t -> SVectorNT{3O, N}( _wrap_vector6(cache, eph, parentid, id, t) ),
        t -> SVectorNT{3O, N}( _wrap_vector9(cache, eph, parentid, id, t) ),
        t -> SVectorNT{3O, N}( _wrap_vector12(cache, eph, parentid, id, t) )
    )

    add_point!(frames, name, id, axesid, POINT_CLASSID_DYNAMIC, funs, parentid) 
    nothing 
 
end

function _wrap_vector3(cache, eph, parentid, id, t)
    tid = Threads.threadid()
    c = cache[tid]
    ephem_compute!(c, eph, 2451545, t / 86400,  id, parentid, 0)
    @inbounds @views return SA[c[1:3]...]
end

function _wrap_vector6(cache, eph, parentid, id, t)
    tid = Threads.threadid()
    c = cache[tid]
    ephem_compute!(c, eph, 2451545, t / 86400,  id, parentid, 1)
    @inbounds @views return SA[c[1:6]...]
end

function _wrap_vector9(cache, eph, parentid, id, t)
    tid = Threads.threadid()
    c = cache[tid]
    ephem_compute!(c, eph, 2451545, t / 86400,  id, parentid, 2)
    @inbounds @views return SA[c[1:9]...]
end

function _wrap_vector12(cache, eph, parentid, id, t)
    tid = Threads.threadid()
    c = cache[tid]
    ephem_compute!(c, eph, 2451545, t / 86400,  id, parentid, 3)
    @inbounds @views return SA[c[1:12]...]
end

end