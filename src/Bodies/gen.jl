export generate_body!, parse_naifnames

import Basic.Utils: genf_psngin, format_camelcase

const DEF_NAIF2NAME = joinpath(@__DIR__, "..", "..", "res", "naif2name.txt")

"""
    parse_naifnames(file)::Dict{Integer, Symbol}

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
function parse_naifnames(file::String)
    matched = match.(r"(?<id>[0-9]{1,})\W{1,}(?<name>\w{1,})\W{1,}(?<type>\w{1,})", readlines(file))
    return Dict(
        typeof(m) == RegexMatch ? 
            parse(Int, m[:id]) => format_camelcase(Symbol, lowercase(m[:name])) : () 
                for m in matched
    ), Dict(
        typeof(m) == RegexMatch ? 
            parse(Int, m[:id]) => format_camelcase(Symbol, lowercase(m[:type])) : () 
                for m in matched
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
    generate_body!(gen::String, bname::Symbol, bid::N, 
        objtype::Symbol, data::D) where {N<:Integer, D<:AbstractDict}

Generate properties for a celestial object. 

### Input/s 

- `gen` -- In-place modified generated code 
- `bid` -- Body NAIF id 
- `objtype` -- Body object type. Shall be one of `AbstractBody` subtypes.
- `data` -- Bodies data dictionary of type `Dict{NAIFId, Any}`

### Output 

`gen` is modified and given as output.
"""
function generate_body!(gen::String, bname::Symbol, bid::N, 
    objtype::Symbol, data::D) where {N<:Integer, D<:AbstractDict}
    
    @debug "[Bodies] Generating code for $bname"

    sid = "$(bid)"
    if objtype == :NaturalSatellite || objtype == :Planet
        bseb = sid[1]
    elseif objtype == :Barycenter 
        bseb = bid 
    else 
        bseb = 0
    end

    # parse code 
    gen *= "\n"
    gen *= "#%BODIES::$bname\n"
    # subtype 
    gen *= "struct $(bname) <: $(objtype) end\n"
    gen *= genf_psngin(:Bodies, :body_naifid, bid, (nothing, "Type{$bname}"))
    gen *= genf_psngin(:Bodies, :body_from_naifid, bname, (nothing, Val{bid}))
    if bid != 0
        gen *= genf_psngin(:Bodies, :body_gm, data[bid][:gm], (nothing, Val{bid}))
    else 
        gen *= genf_psngin(:Bodies, :body_gm, 
            sum([data[i][:gm] for i in range(1,10)]), (nothing, Val{bid}))
    end
    gen *= genf_psngin(:Bodies, :body_system_equivalent, bseb, (nothing, Val{bid}))
    if bid > 9 
        gen *= genf_psngin(:Bodies, :body_equatorial_radius, 
            data[bid][:radii][1], (nothing, Val{bid}))
        gen *= genf_psngin(:Bodies, :body_polar_radius, 
            data[bid][:radii][3], (nothing, Val{bid}))
        gen *= genf_psngin(:Bodies, :body_mean_radius, 
            sum(data[bid][:radii])/length(data[bid][:radii]), 
            (nothing, Val{bid}))
    else # barycenters 
        prop = (
            :body_equatorial_radius,
            :body_polar_radius, :body_mean_radius
        )
        for p in prop 
            sp = join(uppercasefirst.(split.("$(p)", "_")[2:end]), " ")
            gen *= genf_psngin(:Bodies, p, 
                """throw(error("[Bodies] $(bname) has no $sp"))""", 
                (nothing, Val{bid}))
        end
    end
    gen
end
