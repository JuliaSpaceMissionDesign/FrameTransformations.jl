const DEF_NAIF2NAME = joinpath(@__DIR__, "..", "..", "res", "naif2name.txt")

"""
    naif2names()::Dict{Integer, Symbol}

Parse NAIFId-object name equivalence dictionary from `file`.

The supplied `file` must be written in the following format:
```text
0           'SOLAR_SYSTEM_BARYCENTER'
499         'MARS'
```
i.e. first the NAIFId, than the body name using the SPICE convention, i.e. 
uppercase, snakecase names. The name must be included between quotation marks, 
as in the example.
"""
function naif2names(file::String)
    matched = match.(r"(?<id>[0-9]{1,})\W{1,}(?<name>\w{1,})", readlines(file))
    return Dict(
        typeof(m) == RegexMatch ? 
            parse(Int, m[:id]) => format_camelcase(Symbol, lowercase(m[:name])) : () 
                for m in matched
    )
end
""" 
    naif2names()::Dict{Integer, Symbol}

Parse NAIFId-object name equivalence dictionary from `res` folder default file.
"""
function naif2names()
    naif2names(DEF_NAIF2NAME)
end