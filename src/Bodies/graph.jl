const BODIES = NodeGraph{NAIFId, NAIFId}(SimpleGraph())

register_body!(id) = add_vertex!(BODIES, id)
connect_bodies!(id1, id2) = add_edge!(BODIES, id1, id2)
find_path(from::CelestialBody, to::CelestialBody) = get_edgenodes(BODIES, body_naifid(from), body_naifid(to))

const BODY_TYPES = (
    :CelestialBody,
    :Barycenter,
    :Planet,
    :NaturalSatellite,
    :MinorBody,
    :Asteroid,
    :Comet,
)

macro add_body(name::Symbol, id::Int, super::Symbol, args...)
    name_str = String(name)
    typ_str = join(uppercasefirst.(split(name_str, '_')))
    typ = Symbol(typ_str)
    parent = nothing
    barycenter = nothing
    _export = false
    if !(super in BODY_TYPES)
        throw(ArgumentError("Invalid supertype: $super"))
    end
    for a in args
        a isa Expr || continue
        a.head == :(=) || continue

        if a.args[1] == :name
            val = a.args[2]
            val isa Symbol || throw(ArgumentError("Invalid argument: $a"))
            typ = val
            typ_str = String(typ)
        elseif a.args[1] == :parent
            val = a.args[2]
            val isa Symbol || throw(ArgumentError("Invalid argument: $a"))
            parent = val
        elseif a.args[1] == :barycenter
            val = a.args[2]
            val isa Symbol || throw(ArgumentError("Invalid argument: $a"))
            barycenter = val
        elseif a.args[1] == :do_export
            val = a.args[2]
            val isa Bool || throw(ArgumentError("Invalid argument: $a"))
            _export = val
        end
    end
    parts = split(String(typ), r"(?=[A-Z1-9])")
    show_str = join(parts, " ")
    doc_str = parts[end] in ("Barycenter", "Sun") ? join(["the"; parts], " ") : show_str
    super_str = String(super)
    reg = if parent === nothing
        :(Bodies.register_body!($id))
    else
        :(Bodies.connect_bodies!(body_naifid($parent), $id))
    end
    exp = _export ? :(export $name, $typ) : :()
    parent_expr = parent !== nothing ? :(body_parent(::$typ) = $parent) : :()
    bseb_expr = :(body_system_equivalent(::$typ) = $barycenter)
    id_expr = :(Bodies.body_naifid(::$typ) = $id)
    fromid_expr = :(Bodies.from_naifid(::Val{$id}) = $name)
    return quote
        """
            $($typ_str) <: $($super_str)
        A type representing $($doc_str).
        """
        struct $(esc(typ)) <: $(esc(super)) end

        """
            $($name_str)
        The singleton instance of the [`$($typ_str)`](@ref) type.
        """
        const $(esc(name)) = $(esc(typ))()
        Base.show(io::IO, ::$(esc(typ))) = print(io, "$($show_str)")
        $(esc(id_expr))
        $(esc(fromid_expr))
        $(esc(reg))
        $(esc(parent_expr))
        $(esc(bseb_expr))
        $(esc(exp))
        nothing
    end
end

@add_body ssb 0     Barycenter      name=SolarSystemBarycenter  do_export=true
@add_body sun 10    CelestialBody   parent=ssb                  do_export=true
