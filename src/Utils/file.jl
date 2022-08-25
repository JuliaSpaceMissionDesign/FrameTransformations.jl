abstract type AbstractFile end

const FILEFORMAT = (
    :TPC, 
    :JSON,
)

for fmt in FILEFORMAT
    @eval begin
        @make_struct_fromschema $(fmt) $AbstractFile (path, String)
        export $(fmt)
    end
end

filepath(file::T) where {T<:AbstractFile} = file.path