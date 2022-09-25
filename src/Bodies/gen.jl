export generate_body!

import Basic.Utils: genf_psngin, format_camelcase


"""
    generate_body!(gen::String, bname::Symbol, bid::N, 
        objtype::Symbol, data::D) where {N<:Integer, D<:AbstractDict}

Generate properties for a celestial object. 

### Input/s 

- `gen` -- In-place modified generated code 
- `bid` -- Body NAIF id 
- `bname` -- Body name
- `objtype` -- Body object type. Shall be one of `AbstractBody` subtypes.
- `data` -- Bodies data dictionary containing the properties to parse

### Output 

`gen` is modified and given as output.
"""
function generate_body!(gen::String, bname::Symbol, bid::N, 
    objtype::Symbol, data::D) where {N<:Integer, D<:AbstractDict}

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
    gen *= """
           struct $(bname) <: $(objtype) end
           """
    gen *= genf_psngin(:Bodies, :body_naifid, "NAIFId($bid)", (nothing, bname))
    # gen *= genf_psngin(:Bodies, :body_from_naifid, bname, (nothing, Val{bid})) # removed 
    if haskey(data, :gm)
        gen *= genf_psngin(:Bodies, :body_gm, pop!(data, :gm), (nothing, Val{bid}))
    else 
        throw(error("[Bodies] gravitational parameter must be loaded for $bname to declare the body!"))
    end
    gen *= genf_psngin(:Bodies, :body_system_equivalent, "NAIFId($bseb)", (nothing, Val{bid}))
    if bid > 9 
        for (k, v) in data
            gen *= genf_psngin(:Bodies, Symbol("body_", k), v, (nothing, Val{bid}))
        end
    else # barycenters 
        for (k, _) in data
            sp = join(uppercasefirst.(split.("$(k)", "_")[2:end]), " ")
            gen *= genf_psngin(:Bodies, Symbol("body_", k), 
                """throw(error("[Bodies] $(bname) has no $sp"))""", 
                (nothing, Val{bid}))
        end
    end
    gen
end
