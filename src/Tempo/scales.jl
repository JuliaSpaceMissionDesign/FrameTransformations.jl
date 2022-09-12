using NodeGraphs
import Basic: register!, connect!, find_path

export TIMESCALES

struct NotATimeScale <: TimeScale end

const NAMES = (
    # :UniversalTime,
    :InternationalAtomicTime,
    :TerrestrialTime,
    :GeocentricCoordinateTime,
    :BarycentricCoordinateTime,
    :BarycentricDynamicalTime,
)

const ACRONYMS = (
    # :UT1,
    :TAI,
    :TT,
    :TCG,
    :TCB,
    :TDB,
)

for (acronym, name) in zip(ACRONYMS, NAMES)
    acro_str = String(acronym)
    name_str = String(name)
    name_split = join(split(name_str, r"(?=[A-Z])"), " ")
    wiki = replace(name_split, " "=>"_")
    @eval begin
        """
            $($name_str)
        A type representing the $($name_split) ($($acro_str)) time scale.
        # References
        - [Wikipedia](https://en.wikipedia.org/wiki/$($wiki))
        """
        struct $name <: TimeScale end

        """
            $($acro_str)
        The singleton instance of the [`$($name_str)`](@ref) type representing
        the $($name_split) ($($acro_str)) time scale.
        # References
        - [Wikipedia](https://en.wikipedia.org/wiki/$($wiki))
        """
        const $acronym = $name()

        export $name, $acronym

        Base.show(io::IO, ::$name) = print(io, "$($acro_str)")
        tryparse(::Val{Symbol($acro_str)}) = $acronym
    end
end


"""
    TIMESCALES = NodeGraph{TimeScale, UInt8, Int16}(SimpleGraph{UInt8}())

Time Scales Graph singleton.
"""
const TIMESCALES = NodeGraph{TimeScale, UInt8, Int16}(SimpleGraph{UInt8}())

"""
    register!(s::TS) where {TS<:TimeScale}

Register a new `TimeScale` subtype to the graph.
"""
function register!(s::TS) where {TS<:TimeScale}
    add_vertex!(TIMESCALES, s)
end

"""
    connect!(s1::TS1, s2::TS2; oneway=false) where {TS1<:TimeScale, TS2<:TimeScale}

Connect two new `TimeScale` subtype in the graph.
"""
function connect!(s1::TS1, s2::TS2; 
    oneway=false) where {TS1<:TimeScale, TS2<:TimeScale} 
    add_edge!(TIMESCALES, s1, s2)
    oneway || add_edge!(TIMESCALES, s2, s1)
end

connect!(TAI, TT)
# connect!(TAI, UT1)
connect!(TT, TCG)
connect!(TT, TDB)
connect!(TCB, TDB)

function find_path(from::TS1, to::TS2) where {TS1<:TimeScale, TS2<:TimeScale} 
    get_nodes(TIMESCALES, from, to)
end