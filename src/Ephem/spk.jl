using Calceph: Ephem as CalcephEphem, unsafe_compute!, useNaifId, unitKM, unitSec

struct SPK <: AbstractEphemeris
    handler::CalcephEphem
    function SPK(files::Vector{<:AbstractString})
        new(CalcephEphem(files))
    end
end

SPK(file::AbstractString) = SPK([file])

# TODO: check if this offset is right
to_julian(jd2000::Float64) = jd2000 + 2451545.0

function vector3!(r::Vector{Float64}, eph::SPK, from::Int, to::Int, jd2000::Float64; time::Float64=0.0)
    unsafe_compute!(r, eph.handler, to_julian(jd2000), time, from, to, useNaifId+unitKM+unitSec, 0)
end

function vector6!(rv::Vector{Float64}, eph::SPK, from::Int, to::Int, jd2000::Float64; time::Float64=0.0)
    unsafe_compute!(rv, eph.handler, to_julian(jd2000), time, from, to, useNaifId+unitKM+unitSec, 1)
end

function vector9!(rva::Vector{Float64}, eph::SPK, from::Int, to::Int, jd2000::Float64; time::Float64=0.0)
    unsafe_compute!(rva, eph.handler, to_julian(jd2000), time, from, to, useNaifId+unitKM+unitSec, 2)
end

function vector12!(rvaj::Vector{Float64}, eph::SPK, from::Int, to::Int, jd2000::Float64; time::Float64=0.0)
    unsafe_compute!(rvaj, eph.handler, to_julian(jd2000), time, from, to, useNaifId+unitKM+unitSec, 3)
end