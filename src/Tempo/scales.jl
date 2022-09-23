using NodeGraphs
import Basic: register!, connect!, find_path

export TIMESCALES

struct NotATimeScale <: TimeScale end

const TIMESCALES_NAMES = (
    # :UniversalTime,
    :InternationalAtomicTime,
    :CoordinatedUniversalTime,
    :TerrestrialTime,
    :GeocentricCoordinateTime,
    :BarycentricCoordinateTime,
    :BarycentricDynamicalTime,
)

const TIMESCALES_ACRONYMS = (
    # :UT1,
    :TAI,
    :UTC,
    :TT,
    :TCG,
    :TCB,
    :TDB,
)

for (acronym, name) in zip(TIMESCALES_ACRONYMS, TIMESCALES_NAMES)
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
        Base.tryparse(::Val{Symbol($acro_str)}) = $acronym
    end
end

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
connect!(TAI, UTC)
# connect!(TAI, UT1)
connect!(TT, TCG)
connect!(TT, TDB)
connect!(TCB, TDB)

function find_path(from::TS1, to::TS2) where {TS1<:TimeScale, TS2<:TimeScale} 
    get_nodes(TIMESCALES, from, to)
end


"""
    TIMESCALES

Time Scales Graph singleton. 

It is a `NodeGraph{TimeScale, UInt8, Int16}` type instance, containing the 
following scales: $(String.(TIMESCALES_ACRONYMS))

It can be easily extended using the `register!` and/or `connect!` methods and 
defining an appropriate type for the new item, as well as its relation with 
the other nodes in the graph.

### Example

```@example
# define a new timescale type - shall be a TimeScale subtype
struct NewTimeScale <: TimeScale end 

# create the user-defined timescale singleton instance
const NTS = NewTimeScale()

# define offset to and from another timescale in the graph 
offset(::InternationalAtomicTime, ::NewTimeScale) = 1.0
offset(::NewTimeScale, ::InternationalAtomicTime) = -1.0

# connect to the graph, with the parent node used in `offset`
connect!(TAI, NTS)
```
"""
TIMESCALES