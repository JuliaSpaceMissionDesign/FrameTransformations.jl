export parse_universe

using Basic.Utils: YAML, JSON, filepath
import YAML as YAMLLib
import JSON3

function _parse_universe(configfile::YAML)
    YAMLLib.load_file(filepath(configfile); dicttype=OrderedDict{String, Any})
end

function _parse_universe(configfile::JSON) 
    open(filepath(configfile), "r") do f 
        data = JSON3.read(f, OrderedDict{String, Any})
        return data
    end
end

function parse_universe(configfile::Union{YAML, JSON})
    config = _parse_universe(configfile)
    data = OrderedDict()
    push!(data, :gen => Vector{Tuple{GenMeta, String}}())

    # validation step 
    if isvaliduniverse(config)
        for k in keys(config)
            parser = Symbol("parse_$(k)!")
            @info "[Universe] Parsing $k..."
            @eval begin 
                $parser($data, $config)
            end
        end
    else
        # if not valid, throw error
        validate(UniverseSchema, config)
        throw(error("[Universe] Input universe file is not valid!"))
    end
    return data 
end