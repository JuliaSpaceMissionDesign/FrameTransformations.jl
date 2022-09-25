export parse_bodies!

using Basic.Bodies: generate_body!, NAIFId


function parse_bodies!(datadict::D1, 
    configdict::D2) where {D1 <: AbstractDict, D2 <: AbstractDict}
    field = configdict["bodies"]
    bgen = ""

    bdata = OrderedDict{NAIFId, Symbol}()
    for body in field
        name, bid, bgen = _parse_body!(bgen, datadict, body)
        push!(
            bdata, 
            NAIFId(bid) => name
        )
    end
    push!(
        datadict,
        :bodies => OrderedDict(
            :naifid2name => bdata, 
            :name2naifid => OrderedDict(v => k for (k, v) in bdata)
        )
    )
    bcodeid = bytes2hex(sha256(bgen))
    push!(
        datadict[:gen], 
        (GenMeta("Universe/Bodies", bcodeid), bgen)
    )
    @info "[Universe] Autogen code for bodies with id: $bcodeid"
    nothing
end

function _parse_body!(gen::String, data::D1, 
    args::D2) where {D1<:AbstractDict, D2<:AbstractDict}

    name, point = format_camelcase(Symbol, args["name"]), args["point"]

    # body id, name and type 
    bid, bname, bobtype = __get_body_identifiers(data[:mappings][:naifnames], point)

    # body constants
    bconst = __get_body_constants(data[:constants], bid, args)

    # body models
    bmodels = __get_bodymodeldescriptor(args)

    @info "[Universe/Bodies] Generating code for $bname ($bid)..."
    # generate body 
    gen = generate_body!(gen, name, bid, bobtype, bconst)

    # generate model associated to body
    if bmodels !== nothing
        for (mod, prop) in bmodels 
            for (k, v) in prop    
                gen *= genf_sngin(Symbol("body_", mod, "_", k), 
                    k == :model ? ":$v" : v, (nothing, name))
            end
        end
    end
    return name, bid, gen
end

function __get_body_identifiers(data::D, bname::String) where {D<:AbstractDict}
    bname = Symbol(bname)
    bid = data[:name2id][bname]
    bobtype = data[:id2type][bid]
    return bid, bname, bobtype
end

function __get_body_identifiers(data::D, bid::Int) where {D<:AbstractDict}
    bname = data[:id2name][bid]
    bobtype = data[:id2type][bid]
    return bid, bname, bobtype
end

function __get_body_constants(data::D1, bid::Int, 
    args::D2) where {D1<:AbstractDict, D2<:AbstractDict}

    constants = OrderedDict{Symbol, Float64}()
    push!(constants, :gm => haskey(args, "gm") ? args["gm"] : data[bid][:gm])
    if haskey(args, "radius")
        eqrad = args["radius"]
        f = haskey(args, "flattening") ? args["flattening"] : 0.0
        polrad = eqrad/(1+f)
        meanrad = (eqrad + polrad)/2
    else 
        if haskey(data[bid], :radii)
            eqrad = data[bid][:radii][1]
            polrad = data[bid][:radii][3]
            meanrad = sum(data[bid][:radii])/3
            f = (eqrad - polrad)/polrad
        else
            eqrad, polrad, meanrad, f = NaN, NaN, NaN, NaN
        end
    end
    push!(
        constants, 
        :equatorial_radius => eqrad,
        :polar_radius => polrad, 
        :mean_radius => meanrad,
        :flattening => f
    )
    return constants
end

function __get_bodymodeldescriptor(args)
    if haskey(args, "models")
        # handle models 
        modelmap = OrderedDict{Symbol, Dict{Symbol, String}}()
        for model in args["models"]
            mtype = Symbol(pop!(model, "typed"))
            push!(
                modelmap,
                mtype => OrderedDict{Symbol, String}(Symbol(k) => v for (k, v) in model)
            )
        end
        return modelmap
    end 
    nothing 
end
