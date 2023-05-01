module Tempo

import FunctionWrappers: FunctionWrapper
using RemoteFiles: @RemoteFile, download
# TODO: remove this dependency - could be handled by Tempo itself?
using Dates: DateTime as DatesDateTime, datetime2julian, now

using Basic

using InterfacesUtils.Utils: format_camelcase
using InterfacesUtils.Math: interpolate
using InterfacesUtils.Interfaces.Errors: AbstractGenericException, @module_error

using MultiGraphs:
    MappedNodeGraph,
    MappedDiGraph,
    AbstractGraphNode,
    SimpleDiGraph,
    has_vertex,
    add_edge!,
    get_path,
    get_mappedid,
    get_mappednode

import MultiGraphs: get_node_id, add_vertex!

import PrecompileTools

const DAY2SEC = 86400.0
const YEAR2SEC = 60.0 * 60.0 * 24.0 * 365.25
const CENTURY2SEC = 60.0 * 60.0 * 24.0 * 365.25 * 100.0
const CENTURY2DAY = 36525.0

"""
    DJ2000

Reference epoch (J2000.0), Julian Date (`2451545.0`). 
It is `12:00 01-01-2000`.
"""
const DJ2000 = 2451545.0

"""
    DMJD

Reference epoch (J2000.0), Modified Julian Date (`51544.5`).
"""
const DMJD = 51544.5

"""
    DJM0

Julian Date of Modified Julian Date zero point (`2400000.5`).
It is `00:00 17-11-1858`.
"""
const DJM0 = 2400000.5

"""
AbstractTimeScale

All timescales are subtypes of the abstract type `AbstractTimeScale`.
"""
abstract type AbstractTimeScale end

export DJ2000, DMJD, DJM0

include("errors.jl")
include("convert.jl")
include("parse.jl")

include("leapseconds.jl")
include("offset.jl")
include("scales.jl")

include("datetime.jl")
include("origin.jl")
include("epoch.jl")

# Precompilation routines 
PrecompileTools.@setup_workload begin 

    PrecompileTools.@compile_workload begin 

        # Precompile Epochs routines for all time scales 
        for scale in Tempo.TIMESCALES_ACRONYMS
            epo = Epoch("2022-02-12T12:00:34.3241 $scale")
            j2000(epo)
            j2000s(epo)
            j2000c(epo)
        end 

        # Precompile timescale offset functions 
        for fcn in (offset_gps2tai, offset_tai2gps, offset_tt2tdbh, offset_utc2tai, 
                    offset_tai2utc, offset_tdb2tt, offset_tt2tdb, offset_tdb2tcb, 
                    offset_tcb2tdb, offset_tt2tcg, offset_tcg2tt, offset_tt2tai, 
                    offset_tai2tt)

            fcn(200.0)
        end

        # Precompile datetime stuff 
        Date(23)
        DateTime("2022-02-12T12:00:32")
        DateTime(21321.034)
        
        date = Date(2022, 12)
        dt = DateTime(date, 2312.04)

        Date(dt)
        Time(dt)
        
        for fcn in (year, month, day, hour, minute, second)
            fcn(dt)
        end

        # Precompile smaller routines
        tai2utc(DJ2000, 0.0)
        utc2tai(DJ2000, 0.0)
        jd2calhms(DJ2000, 0.0)
        fd2hmsf(0.4)
    end
end


end
