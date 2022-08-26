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
    nothing
end

function write(file::JSON, data) 
    write_json(filepath(file), data)
end

function write(file::ADF, data::AbstractDict{String, Any}; 
    mode=:write, meta=nothing)
    jldopen(filepath(file), FILE_OPENMODE[mode]) do file 
        for (k, v) in pairs(data)
            file[k] = v
        end
        meta !== nothing ? file["meta"] = meta : nothing
    end
    nothing
end