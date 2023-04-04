export write

import Base: write

"""
    write(filename::T, data::Any; kwargs...) where T <: AbstractFile

Generic writer for different file/s format.
"""
function write() end

function write_json(path::String, data)
    open(path, "w") do f
        JSON3.pretty(f, data)
    end
    return nothing
end

function write(file::JSON, data)
    return write_json(filepath(file), data)
end
