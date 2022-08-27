using Basic.Utils: format_camelcase

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

template_inlineconst(body, fun, value) = "@inline $fun(::$body) = $value\n"
template_structwithparent(name, parent) = "struct $(name) <: $(parent) end\n"
template_singleton(name, typ) = "const $name = $typ()\n"

function template_inlineconst(typs::Tuple, fun, value)
    string = ""
    for typ in typs
        string *= template_inlineconst(typ, fun, value)
    end
    string 
end

"""
    generate_body!(gen::String, gname::Symbol, bname::Symbol, 
        bid::N, centid::Union{N, Nothing}, objtype::Symbol, 
        data::D) where {N<:Integer, D<:AbstractDict}

Generate properties for a celestial object. 

### Input/s 

- `gen` -- In-place modified generated code 
- `gname` -- Name of the graph where to append the objects 
- `bname` -- Body NAIF name 
- `bid` -- Body NAIF id 
- `centid` -- Body center id 
- `objtype` -- Body object type. Shall be one of `AbstractBody` subtypes.
- `data` -- Bodies data dictionary of type `Dict{NAIFId, Any}`

### Output 

`gen` is modified and given as output.
"""
function generate_body!(gen::String, gname::Symbol, 
    bname::Symbol, bid::N, centid::Union{N, Nothing}, 
    objtype::Symbol, data::D) where {N<:Integer, D<:AbstractDict}

    sid = "$(bid)"
    if objtype == :NaturalSatellite || objtype == :Planet
        bseb = sid[1]
    elseif objtype == :Barycenter 
        bseb = nothing 
    else 
        bseb = 0
    end
    bnameid = (bname, Val{bid})

    # parse code 
    gen *= "\n"
    gen *= "#%$bname\n"
    gen *= template_structwithparent(bname, objtype)
    if centid === nothing
        gen *= "register!($gname, $bid)\n"
    else 
        gen *= "connect!($gname, $centid, $bid)\n"
    end
    gen *= template_inlineconst(bname, :body_naifid, bid)
    gen *= template_inlineconst(bnameid, :body_parent, 
        centid === nothing ? 0 : centid)
    # gen *= template_inlineconst(Val{bid}, :body_from_naifid, bname)
    if bid != 0
        gen *= template_inlineconst(bnameid, :body_gm, data[bid][:gm])
    else 
        gen *= template_inlineconst(bnameid, :body_gm, 
            sum([data[i][:gm] for i in range(1,10)]))
    end
    if bid > 9 
        gen *= template_inlineconst(bnameid, :body_system_equivalent, bseb)
        gen *= template_inlineconst(bnameid, :body_equatorial_radius, 
            data[bid][:radii][1])
        gen *= template_inlineconst(bnameid, :body_polar_radius, 
            data[bid][:radii][3])
        gen *= template_inlineconst(bnameid, :body_mean_radius, 
            sum(data[bid][:radii])/length(data[bid][:radii]))
    else # barycenters 
        prop = (
            :body_system_equivalent, :body_equatorial_radius,
            :body_polar_radius, :body_mean_radius
        )
        for p in prop 
            sp = join(uppercasefirst.(split.("$(p)", "_")[2:end]), " ")
            gen *= template_inlineconst(bname, 
                p, 
                """throw(error("[Bodies] $(bname) has no $sp"))""")
            gen *= template_inlineconst(Val{bid}, 
                p, 
                """throw(error("[Bodies] $(bname) has no $sp"))""")
        end
    end
    gen
end
