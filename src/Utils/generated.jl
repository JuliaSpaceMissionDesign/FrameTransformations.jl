export GEN, GenMeta, parse_generated, template_generated, template_metadata

"""
    GenMeta 

Generated file metadata.

### Fields

- `generated` -- module from which the file is generated 
- `id` -- unique id of the generated body
"""
struct GenMeta
    generated::String 
    id::String
end

Base.isequal(m1::GenMeta, m2::GenMeta) = m1.id == m2.id && m1.generated == m2.generated

"""
    GEN <: AbstractFile

Generated file.

### Fields 

- `path` -- path of the file 
- `meta` -- is a vector of `GenMeta` types 
- `body` -- vector of the bodies of the code corresponding to `meta`s.

### Constructors 

- `GEN(path, meta, body)`
- `GEN(path)`
"""
struct GEN <: AbstractFile 
    path::String
    meta::Vector{GenMeta}
    body::Vector{String}
end
function GEN(path::String)
    m, d = parse_generated(read_generated(path))
    GEN(path, m, d)
end

"""
    load(file::GEN)::GEN
Open a GEN file and parse its data.
"""
function load(file::GEN)
    return file
end

read_generated(file::String) = join(load(TEXT(file)), "\n")

"""
    parse_generated(data::String)

Parse a string with `astronaut` generated file format into a `GEN` type.
"""
function parse_generated(data::String)
    # find all metadata groups
    metadata = (data[id] 
        for id in findall(r"(?<=#%META\n)[\S\s]+?(?=#%ENDS)", data))
    meta = Vector{GenMeta}()
    # parse metadata from the generated file 
    d = Dict{Symbol, String}()
    for mt in metadata 
        for m in split(mt, "\n")[1:end-1]
            sm = split(m, " ")
            push!(d, Symbol(sm[2]) => sm[3])
        end
        push!(meta, GenMeta(d[:generated], d[:id]))
    end

    # find all code groups 
    body = [data[id] 
        for id in findall(r"(?<=#%BODY\n)[\S\s]+?(?=#%ENDS)", data)]
    return meta, body
end

"""
    template_metadata(by::String, id::String) 

Template to parse generated files metadata
"""
function template_metadata(by::String, id::String) 
    """
    #%META
    # generated $(by)
    # id $(id)
    #%ENDS
    """
end

"""
    template_generated(code::String) 

Template to parse generated code
"""
function template_generated(code::String)
    """
    #%BODY
    $(code)
    #%ENDS

    """
end