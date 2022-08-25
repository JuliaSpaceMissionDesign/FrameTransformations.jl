export write

import Base: write

"""
    load(filename::T, data::Any) where T <: AbstractFile

Generic writer for different file/s format.
"""
function write() end 

function write(file::JSON, data) 
    write_json(filepath(file), data)
end

function write_json(path::String, data)
    open(path, "w") do f
        JSON3.pretty(f, data)
    end
    nothing
end