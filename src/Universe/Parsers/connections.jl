export parse_connections!

function parse_connections!(datadict::D1, 
    configdict::D2) where {D1 <: AbstractDict, D2 <: AbstractDict}
    field = configdict["connections"]
    gen = ""
    for group in field
        gen = _parse_connections_graph(gen::String, group["name"], 
            group["points"], datadict)
    end
    codeid = bytes2hex(sha256(gen))
    push!(
        datadict[:gen], 
        (GenMeta("Universe/Connections", codeid), gen)
    )
    @info "[Universe/Connections] Autogen code for `connections` with id: $codeid"
end

function _parse_connections_graph(gen::String, name, points, data)

    @info "[Universe/Connections] parsing $name graph of points..."
    gname = Symbol(format_camelcase(Symbol, String(name)), "Graph")
    gen *= "#%CONNECTIONS::$name\n"
    gen *= "const $gname = Bodies.BodyGraph(Int64, CelestialBody)\n"
    for point in points
        this = point["name"]
        if !haskey(data[:bodies], Symbol(this))
            throw(error(""))
        end
        parent = point["parent"]
        if strip(parent) == ""
            # ssb 
            parent = :SolarSystemB
        end
        gen *= "@inline Bodies.body_parent(::$(Symbol(this, "Type"))) = $parent\n"
        gen *= "connect!(gname, $this, $parent)\n"
    end
    return gen
end