export parse_universe

using Basic.Utils: YAML, JSON, filepath
import YAML as YAMLLib
import JSON3

const DEF_UNIVERSE_FILE = joinpath(@__DIR__, "..", "..", "gen/src", "temp.jl")
const UNIVERSE_EVAL_SEQ = (
    :mappings,
    :constants, 
    :ephemeris, 
    :bodies, 
    :connections, 
    :frames
)

function _parse_universe(configfile::YAML)
    YAMLLib.load_file(filepath(configfile); dicttype=OrderedDict{String, Any})
end

function _parse_universe(configfile::JSON) 
    open(filepath(configfile), "r") do f 
        data = JSON3.read(f, OrderedDict{String, Any})
        return data
    end
end

function parse_universe(configfile::Union{YAML, JSON}, file=DEF_UNIVERSE_FILE)
    config = _parse_universe(configfile)
    # insert generated configfile at the top 
    raw = join(["#@config $line" for line in Base.readlines(filepath(configfile))], "\n")
    fileid = bytes2hex(sha256(raw))

    data = OrderedDict()
    push!(data, :gen => Vector{Tuple{GenMeta, String}}())
    push!(data[:gen], (GenMeta("Universe/@config", fileid), raw))

    # validation step 
    if isvaliduniverse(config)
        for k in UNIVERSE_EVAL_SEQ
            if haskey(config, String(k))
                parser = Symbol("parse_$(k)!")
                @info "[Universe] Parsing $k..."
                @eval begin 
                    $parser($data, $config)
                end
            end
        end
    else
        # if not valid, throw error
        validate(UniverseSchema, config)
        throw(error("[Universe] Input universe file is not valid!"))
    end
    if length(data[:gen]) > 0
        gen = _parse_generated(file, data[:gen])  
        write(gen)
    end
    return data 
end

function _parse_generated(filename, data)
    GEN(
        filename, 
        [meta for (meta, _) in data],
        [gen for (_, gen) in data],
    )
end

function writelines(gen::GEN)
    s = ""
    for (m, d) in zip(gen.meta, gen.body)
        meta = template_metadata(m.generated, m.id)
        data = template_generated(d)
        s *= meta*"\n"*data 
    end
    return s
end

function write(gen::GEN)
    lines = split(writelines(gen), "\n")
    open(gen.path, "w") do f 
        for line in lines
            Base.write(f, line)
            Base.write(f, "\n")
        end
    end
    nothing
end