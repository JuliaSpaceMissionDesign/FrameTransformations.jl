abstract type AbstractFile end

using SHA: sha256

const FILE_OPENMODE = Dict([
    :write => "w",
    :append => "a+",
    :read => "r"
])

const FILE_FORMAT = (
    :TPC, 
    :JSON,
    :TXT,
    :YAML
)

for fmt in FILE_FORMAT
    @eval begin
        @make_struct_fromschema $(fmt) $AbstractFile (path, String)
        export $(fmt)
    end
end

const YML = YAML
export YML

filepath(file::T) where {T<:AbstractFile} = file.path

"""
    fileid(file::String)::String
    fileid(file::F) where {F <: AbstractFile}

Get a unique string descriptor for a file - based on sha256.
"""
function fileid(file::F) where {F <: AbstractFile}
    fileid(filepath(file))
end

function fileid(file::String)
    open(file, "r") do file 
        return bytes2hex(sha256(file))
    end
end