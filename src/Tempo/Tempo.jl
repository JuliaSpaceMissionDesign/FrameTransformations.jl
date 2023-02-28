module Tempo

    import FunctionWrappers: FunctionWrapper
    using RemoteFiles
    using Dates: DateTime as DatesDateTime, datetime2julian, now

    using Basic
    using Basic.Utils: format_camelcase, interpolate
    using Basic: AstronautGenericException, @create_module_error

    using MultiGraphs: 
        MappedNodeGraph, MappedDiGraph, AbstractGraphNode, SimpleDiGraph, 
        has_vertex, add_vertex!, add_edge!, get_path, get_mappedid, get_mappednode
    import MultiGraphs: get_node_id

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

end