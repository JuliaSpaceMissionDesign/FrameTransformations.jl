export parse_ephemeris!

using Basic.Ephemeris: CalcephProvider, ephem_position_records

function parse_ephemeris!(datadict::D1, 
    configdict::D2) where {D1 <: AbstractDict, D2 <: AbstractDict}

    field = configdict["ephemeris"]
    data = OrderedDict{Symbol, Any}()
    gen = """
    const EPHEM_GRAPH = Bodies.BodyGraph(Int64, NAIFId)
    """
    for line in field
        if lowercase(line["parser"]) == "calceph"
            # parse ephemeris with calceph provider 
            provider = CalcephProvider(line["files"])
            push!(
                data, 
                :calceph => provider
            )
            gen *= "const EPHEM_CALCEPH = CalcephProvider($(line["files"]))\n"
            gen *= "register!(EPHEM_GRAPH, EPHEM_CALCEPH)\n"
        else 
            throw(
                error("[Universe/Ephemeris] $(line["parser"]) is not a recognized ephemeris parser")
            )
        end
        
    end
    push!(
        datadict,
        :ephemeris => data
    )
    codeid = bytes2hex(sha256(gen))
    push!(
        datadict[:gen], 
        (GenMeta("Universe/Ephemeris", codeid), gen)
    )
    @info "[Universe/Ephemeris] Autogen code for `ephemeris` with id: $codeid"
    nothing
end
