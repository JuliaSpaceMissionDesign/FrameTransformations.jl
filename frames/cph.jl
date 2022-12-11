
const libcalceph::String = "libcalceph";

const cph_au::Int = 1 
const cph_km::Int = 2
const cph_day::Int = 4
const cph_sec::Int = 8
const cph_use_naifid::Int = 32
const cph_euler::Int = 64
const cph_nutation::Int = 128

const CalcephNullPointerError::String = "Ephemeris kernels are not properly initialised!";

struct CalcephException 
    msg::String 
end

Base.show(io::IO, c::CalcephException) = println(io, "CalcephException: $(c.msg)")

# Utilities 
macro _check_pointer(ptr, msg)
    return quote 
        ($(esc(ptr)) == C_NULL) && throw(CalcephException($(esc(msg))))
    end
end

macro _check_status(stat, msg)
    return quote 
        ($(esc(stat)) == 0) && throw(CalcephException($(esc(msg))))
    end
end

macro _check_order(order)
    return quote 
        if !(0 <= $(esc(order)) <= 3) 
            throw(CalcephException("Order must be between 0 and 3.")) 
        end
    end
end

mutable struct EphemerisKernels
    ptr::Ptr{Cvoid}
    files::Vector{String}
end

function EphemerisKernels(files::Vector{<:AbstractString})
    files = unique(files) # Avoids duplicates 

    # Loads kernels 
    ptr = ccall((:calceph_open_array, libcalceph), 
                Ptr{Cvoid}, (Int, Ptr{Ptr{UInt8}}), 
                length(files), files)
    
    @_check_pointer ptr "Unable to open ephemeris file(s)!"
    
    # Register object destructor    
    obj = finalizer(_kernelDestructor, EphemerisKernels(ptr, files))
    _prefetch(obj) 
    return obj 
end

EphemerisKernels(file::String) = EphemerisKernels([file])
EphemerisKernels() = EphemerisKernels(C_NULL, [""]);

function Base.show(io::IO, eph::EphemerisKernels)
    if eph.ptr != C_NULL 
        print(io, "$(length(eph.files)) Ephemeris Kernels loaded:\n")
        for f in eph.files 
            println(io, "  "*f)
        end
    else 
        print(io, "No kernels loaded.")
    end

    nothing 
end

function _kernelDestructor(eph::EphemerisKernels)
    eph.ptr == C_NULL && return
    
    ccall((:calceph_close, libcalceph), 
          Cvoid, (Ptr{Cvoid},), eph.ptr)
    eph.ptr = C_NULL 

    nothing 
end

function _prefetch(eph::EphemerisKernels)
    @_check_pointer eph.ptr CalcephNullPointerError
    stat = ccall((:calceph_prefetch, libcalceph), Cint, 
                 (Ptr{Cvoid},), eph.ptr)
    @_check_status stat "Unable to prefetch ephemerides!"
    nothing
end

function is_threadsafe(eph::EphemerisKernels)
    @_check_pointer eph.ptr CalcephNullPointerError
    ccall((:calceph_isthreadsafe, libcalceph), Cint, 
          (Ptr{Cvoid},), eph.ptr) == 1
end


# Time properties 
function get_timescale(eph::EphemerisKernels)
    # 1 = TDB, 2 = TCB
    @_check_pointer eph.ptr CalcephNullPointerError
    stat = ccall((:calceph_gettimescale, libcalceph), Cint, (Ptr{Cvoid},), eph.ptr)

    @_check_status stat "Unable to retrieve ephemeris timescale!" 
    return stat
end

function get_timespan(eph::EphemerisKernels)
    @_check_pointer eph.ptr CalcephNullPointerError 

    ft = Ref{Cdouble}(0)
    lt = Ref{Cdouble}(0)
    cn = Ref{Cint}(0)
    
    stat = ccall((:calceph_gettimespan, libcalceph), 
                 Cint, (Ptr{Cvoid}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cint}),
                 eph.ptr, ft, lt,cn)

    @_check_status stat "Unable to compute ephemeris"
    return ft[], lt[], cn[]
end

# Data inspection
struct PositionRecord 
    target::Int 
    center::Int 
    start_epoch::Float64
    stop_epoch::Float64 
    frame::Int
end

struct OrientRecord 
    target::Int 
    start_epoch::Float64 
    stop_epoch::Float64 
    frame::Int
end

function get_orient_records_count(eph::EphemerisKernels)
    @_check_pointer eph.ptr CalcephNullPointerError

    ccall((:calceph_getorientrecordcount, libcalceph), 
          Cint, (Ptr{Cvoid},), eph.ptr)
end

function get_position_records_count(eph::EphemerisKernels)
    @_check_pointer eph.ptr CalcephNullPointerError

    ccall((:calceph_getpositionrecordcount, libcalceph), 
          Cint, (Ptr{Cvoid},), eph.ptr)
end

function get_position_records(eph::EphemerisKernels, idx::Int)
    @_check_pointer eph.ptr CalcephNullPointerError

    target = Ref{Cint}(0)
    center = Ref{Cint}(0)
    start_epoch = Ref{Cdouble}(0.0)
    stop_epoch = Ref{Cdouble}(0.0)
    frame = Ref{Cint}(0)

    stat = ccall((:calceph_getpositionrecordindex, libcalceph), 
          Cint, (Ptr{Cvoid}, Cint, Ref{Cint}, Ref{Cint}, 
                 Ref{Cdouble}, Ref{Cdouble}, Ref{Cint}), 
                 eph.ptr, idx, target, center, start_epoch, stop_epoch, frame)

    @_check_status stat "Cannot compute ephemeris position records at index $(idx)!"

    PositionRecord(target[], center[], start_epoch[], stop_epoch[], frame[])
end

function get_position_records(eph::EphemerisKernels)
    posv = Vector{PositionRecord}() 
    for i = 1:get_position_records_count(eph)
        push!(posv, get_position_records(eph, i))
    end

    return posv
end

function get_orient_records(eph::EphemerisKernels, idx::Int)
    @_check_pointer eph.ptr CalcephNullPointerError

    target = Ref{Cint}(0)
    start_epoch = Ref{Cdouble}(0.0)
    stop_epoch = Ref{Cdouble}(0.0)
    frame = Ref{Cint}(0)

    stat = ccall((:calceph_getorientrecordindex, libcalceph), 
          Cint, (Ptr{Cvoid}, Cint, Ref{Cint}, Ref{Cdouble}, 
                 Ref{Cdouble}, Ref{Cint}), 
                 eph.ptr, idx, target, start_epoch, stop_epoch, frame)

    @_check_status stat "Cannot compute ephemeris orient records at index $(idx)!"

    OrientRecord(target[], start_epoch[], stop_epoch[], frame[])
end

function get_orient_records(eph::EphemerisKernels)
    posv = Vector{OrientRecord}() 
    for i = 1:get_orient_records_count(eph)
        push!(posv, get_orient_records(eph, i))
    end

    return posv
end


# Ephemeris computation
function unsafe_compute_order!(stv::AbstractVector, eph::EphemerisKernels, 
                              jd0::Number, time::Number, target::Int, center::Int, 
                              unit::Int, order::Int)

    ccall((:calceph_compute_order, libcalceph), Cint, (Ptr{Cvoid}, Cdouble, 
          Cdouble, Cint, Cint, Cint, Cint, Ref{Cdouble}), eph.ptr, jd0, time, 
          target, center, unit, order, stv)
end

function compute_order(eph::EphemerisKernels, jd0::T, time::T, 
                        target::Int, center::Int, unit::Int, 
                        order::Int) where {T <: Number}

    @_check_pointer eph.ptr CalcephNullPointerError
    @_check_order order 

    # T = typeof(promote(jd0, time)[1])
    stv = Vector{T}(undef, 3*(order+1))
    stat = unsafe_compute_order!(stv, eph, jd0, time, target, 
                center, unit, order)

    @_check_status stat "Unable to compute ephemeris"

    return stv 
end

function unsafe_orient_order!(stv::AbstractVector, eph::EphemerisKernels, 
                              jd0::Number, time::Number, target::Int, 
                              unit::Int, order::Int)

    ccall((:calceph_orient_order, libcalceph), Cint, (Ptr{Cvoid}, Cdouble, 
          Cdouble, Cint, Cint, Cint, Ref{Cdouble}), eph.ptr, jd0, time, 
          target, unit, order, stv)
end

function orient_order(eph::EphemerisKernels, jd0::Number, time::Number, 
                      target::Int, unit::Int, order::Int)

    @_check_pointer eph.ptr CalcephNullPointerError
    @_check_order order 

    stv = Vector{T}(undef, 3*(order+1))
    unsafe_orient_order!(stv, eph, jd0, time, target, unit, order)
    
    @_check_status stat "Unable to compute ephemeris"
    return stv 
end