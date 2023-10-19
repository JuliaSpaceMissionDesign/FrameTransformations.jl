export TPC, load

FilesIO.@filetype TPC FilesIO.AbstractFile

"""
    load(file::TPC{1})

Open a JPL ASCII `.tpc` file and parse its data in a dictionary.
"""
function FilesIO.load(file::TPC{1})
    mapped = Dict{Int64,Dict{Symbol,Union{Float64,Int64,Vector{Float64}}}}()
    load_tpc!(mapped, FilesIO.filepath(file))
    return sort(mapped)
end

"""
    load(files::TPC)

Open a group of JPL ASCII `.tpc` files and parse their data in a dictionary.
"""
function FilesIO.load(files::TPC)
    mapped = Dict{Int64,Dict{Symbol,Union{Float64,Int64,Vector{Float64}}}}()
    for file in files
        load_tpc!(mapped, FilesIO.filepath(file))
    end
    return sort(mapped)
end

function load_tpc!(
    dict::Dict{Int64,Dict{Symbol,Union{Float64,Int64,Vector{Float64}}}}, filename::String
)
    # load and strip lines (remove tabs and spaces)
    # extract lines which are within `\begindata` and `\begintext`
    file = join(strip.(readlines(filename)), "\n")
    parsed = findall(r"(?s)(?<=\\begindata).*?(?=\\begintext)", file)
    substrings = join([replace(file[mi], "\n" => " ") for mi in parsed], "  ")
    # extract lines which actually have data using the `BODY**** =` pattern
    # this is useful to filter the header and eventual trailings.
    # This vector contains a list of `BODY******* =` elements which will be 
    # splitted afterwards
    nameidx = findall(r"(BODY\w{1,}\D{1,}=)", substrings)
    # row data are extracted as between square brackets, the `=` 
    # before the brackets is added
    dataidx = findall.(r"=\D{1,}\(([^()]*)\)|(=\D{1,}([0-9]*\.?[0-9]*))", substrings)

    # data are mapped to a dictionary
    for i in eachindex(nameidx)
        if length(nameidx[i]) > 0
            # extract names using idx and split to separate the `=` sign.
            raw_name = substrings[nameidx[i]]

            # parse naif id 
            naif = parse(Int, match(r"(?<=BODY)(.*?)(?=_)", raw_name).match)
            # parse property name 
            prop = Symbol(lowercase(strip(match(r"(?<=\d_)(.*?)(?==)", raw_name).match)))
            data = split(replace(substrings[dataidx[i]], "D" => "E"))

            if data[2] == "("
                # data is a vector
                value = parse.(Float64, data[3:(end - 1)])
            else
                try
                    value = [parse(Int64, data[2])]
                catch
                    value = [parse(Float64, data[2])]
                end
            end
            # create temporary dictionary
            temp = Dict(naif => Dict(prop => length(value) > 1 ? value : value[1]))
            # merge with the input dictionary
            mergewith!(merge!, dict, temp)
        end
    end

    return nothing
    
end
