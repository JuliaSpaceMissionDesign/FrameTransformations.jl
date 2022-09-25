export parse_naifnames, parse_mappings!

const DEF_NAIF2NAME = joinpath(@__DIR__, "..", "..", "res", "naif2name.txt")

"""
parse_naifnames(file)::Dict{Integer, Symbol}

Parse NAIFId-object name equivalence dictionary from `file`.

The supplied `file` must be written in the following format:
```text
0           'SOLAR_SYSTEM_BARYCENTER'       BARYCENTER 
499         'MARS'                          PLANET
```
i.e. first the NAIFId, than the body name using the SPICE convention, i.e. 
uppercase, snakecase names. The name must be included between quotation marks, 
as in the example.
"""
function parse_naifnames(file::String)
    # parse naifnames from file 
    matched = match.(r"(?<id>[0-9]{1,})\W{1,}(?<name>\w{1,})\W{1,}(?<type>\w{1,})", readlines(file))
    id2name = OrderedDict(
        typeof(m) == RegexMatch ? 
            parse(Int, m[:id]) => format_camelcase(Symbol, lowercase(m[:name])) : () 
                for m in matched
    )
    name2id = OrderedDict(v=>k for (k,v) in id2name)

    id2type = OrderedDict(
        typeof(m) == RegexMatch ? 
            parse(Int, m[:id]) => format_camelcase(Symbol, lowercase(m[:type])) : () 
                for m in matched
    )
    return OrderedDict(
        :id2name => id2name, 
        :id2type => id2type,
        :name2id => name2id
    )
end

""" 
parse_naifnames()::Dict{Integer, Symbol}

Parse NAIFId-object name equivalence dictionary from `res` folder default file.
"""
function parse_naifnames()
    parse_naifnames(DEF_NAIF2NAME)
end

"""
    parse_mappings!(datadict::D1, configdict::D2) where {D1 <: AbstractDict, D2 <: AbstractDict}

Parse mappings by `name`.
This automatically call a `parse_<name>` parser to populate the `:mappings` 
entry of the input dictionary.
"""
function parse_mappings!(datadict::D1, 
    configdict::D2) where {D1 <: AbstractDict, D2 <: AbstractDict}
    parsed = OrderedDict{Symbol, Any}()
    for line in configdict["mappings"]
        parser = Symbol("parse_$(line["name"])")
        for file in line["paths"]
            @eval begin
                push!(
                    $parsed,
                    Symbol($(line["name"])) => $parser($(file))
                )
            end
        end
    end
    push!(
        datadict, 
        :mappings => parsed
    )
    nothing
end
