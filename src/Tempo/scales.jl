using NodeGraphs: NodeGraph, add_edge!, add_vertex!, get_nodes

export
    find_path,
    link_scales!,
    register_scale!

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

const SCALES = NodeGraph{TimeScale}()

function register_scale!(s)
    add_vertex!(SCALES, s)
end

function link_scales!(s1, s2; oneway=false)
    add_edge!(SCALES, s1, s2)
    oneway || add_edge!(SCALES, s2, s1)
end

link_scales!(TAI, TT)
# link_scales!(TAI, UT1)
link_scales!(TT, TCG)
link_scales!(TT, TDB)
link_scales!(TCB, TDB)

find_path(from, to) = get_nodes(SCALES, from, to)