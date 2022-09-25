export parse_bodies!

using Basic.Bodies: generate_body!, NAIFId


function parse_bodies!(datadict::D1, 
    configdict::D2) where {D1 <: AbstractDict, D2 <: AbstractDict}
    field = configdict["bodies"]
    push!(
        field, 
        OrderedDict("name" => "SolarSystemB", "point" => 0)
    )
    gen = ""

    bdata = OrderedDict{Symbol, NAIFId}()
    for body in field
        name, bid, gen = _parse_body!(gen, datadict, body)
        push!(
            bdata, 
            name => NAIFId(bid)
        )
    end
    push!(
        datadict,
        :bodies => bdata
    )
    codeid = bytes2hex(sha256(gen))
    push!(
        datadict[:gen], 
        (GenMeta("Universe/Bodies", codeid), gen)
    )
    @info "[Universe/Bodies] Autogen code for `bodies` with id: $codeid"
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
    btname = Symbol(name, "Type")
    gen = generate_body!(gen, btname, bid, bobtype, bconst)
    gen = __remove_duplicated!(gen, name)

    # generate model associated to body
    if bmodels !== nothing
        for (mod, prop) in bmodels 
            for (k, v) in prop    
                gen *= genf_sngin(Symbol("body_", mod, "_", k), 
                    k == :model ? ":$v" : v, (nothing, name))
            end
        end
    end
    # create singleton
    gen *= "const $name = $(btname)()\n"
    return name, bid, gen
end

function __remove_duplicated!(gen::String, bname)
    s = split(gen, "\n")
    uniques = Vector{String}()
    duplicated = false
    for line in s
        if line != ""
            if !(line in uniques)
                push!(
                    uniques, line
                )
            else
                duplicated = true
            end
        else 
            push!(uniques, "")
        end
    end
    if duplicated 
        @warn "[Universe/Bodies] Removing duplicated items for body $bname..."
    end
    return join(uniques, "\n")
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
