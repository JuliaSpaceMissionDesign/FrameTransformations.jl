export generate_body!

import Basic.Utils: format_camelcase
using CodeGen

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
    gen *= generate_fun_single_withmodule(
        :Bodies, :body_naifid, "NAIFId($bid)", (nothing, bname); 
        prefixes=("inline", )
    )
    if haskey(data, :gm)
        gen *= generate_fun_single_withmodule(
            :Bodies, :body_gm, pop!(data, :gm), (nothing, Val{bid}); 
            prefixes=("inline", )
        )
    else 
        throw(error("[Bodies] gravitational parameter must be loaded for $bname to declare the body!"))
    end
    gen *= generate_fun_single_withmodule(
        :Bodies, :body_system_equivalent, "NAIFId($bseb)", (nothing, Val{bid}); 
        prefixes=("inline", )
    )
    if bid > 9 
        for (k, v) in data
            gen *= generate_fun_single_withmodule(
                :Bodies,  Symbol("body_", k), v, (nothing, Val{bid}); 
                prefixes=("inline", )
            )
        end
    else # barycenters 
        for (k, _) in data
            sp = join(uppercasefirst.(split.("$(k)", "_")[2:end]), " ")
            gen *= generate_fun_single_withmodule(
                    :Bodies,  Symbol("body_", k), 
                    """throw(error("[Bodies] $(bname) has no $sp"))""", 
                    (nothing, Val{bid}); 
                    prefixes=("inline", )
                )
        end
    end
    gen
end
