export parse_ephemeris!

using Basic.Ephemeris: CalcephProvider

function parse_ephemeris!(datadict::D1, 
    configdict::D2) where {D1 <: AbstractDict, D2 <: AbstractDict}

    field = configdict["ephemeris"]
    data = OrderedDict{Symbol, Any}()
    for line in field
        if lowercase(line["parser"]) == "calceph"
            # parse ephemeris with calceph provider 
            push!(
                data, 
                :calceph => CalcephProvider(line["files"])
            )
        else 
            throw(
                error("[Universe] $(line["parser"]) is not a recognized ephemeris parser")
            )
        end
        
    end
    push!(
        datadict,
        :ephemeris => data
    )
    nothing
end
