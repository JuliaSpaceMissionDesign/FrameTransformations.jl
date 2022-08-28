export load

"""
    load(filename::T) where T <: AbstractFile
    load(filenames::Vector{T}) where T <: AbstractFile

Generic loader of different file/s format.
"""
function load() end

"""
    load(file::JSON)

Open a JSON file and parse its data in a dictionary.
"""
function load(file::JSON)
    open(filepath(file), "r") do f 
        data = JSON3.read(f)
        return OrderedDict(data)
    end
end

"""
    load(file::TXT)

Open a TEXT file and parse its data in a list of strings.
"""
function load(file::TXT)
    readlines(filepath(file))
end

"""
    load(file::YAML)

Open a YAML file and parse its data in a dictionary.
"""
function load(file::YAML)
    YAMLLib.load_file(filepath(file); dicttype=OrderedDict{Symbol, Any})
end

"""
    load(file::TPC)

Open a JPL ASCII `.tpc` file and parse its data in a dictionary.
"""
function load(file::TPC)
    mapped = Dict{Int64, Dict{Symbol, Union{Float64, Vector{Float64}}}}()
    load_tpc!(mapped, filepath(file))
    sort(mapped)
end

"""
    load(files::Vector{TPC})

Open a group of JPL ASCII `.tpc` files and parse their data in a dictionary.
"""
function load(files::Vector{TPC})
    mapped = Dict{Int64, Dict{Symbol, Union{Float64, Vector{Float64}}}}()
    for file in files
        load_tpc!(mapped, filepath(file))
    end
    sort(mapped)
end

function load_tpc!(dict::Dict{Int64, Dict{Symbol, 
    Union{Float64, Vector{Float64}}}}, filename::String)
    # load and strip lines (remove tabs and spaces)
    # extract lines which are within `\begindata` and `\begintext`
    parsed = split(join(strip.(readlines(filename)), " "), 
        r"(?<=\\begintext).*?(?=\\begindata\s*BODY*)") 
    # extract lines which actually have data using the `BODY**** =` pattern
    # this is useful to filter the header and eventual trailings.
    # This vector contains a list of `BODY******* =` elements which will be 
    # splitted afterwards
    names_idx = findall.(r"(BODY\w{1,}\D{1,}=)", parsed)
    # row data are extracted as between square brackets, the `=` 
    # before the brackets is added
    datas_idx = findall.(r"=\D{1,}\(([^()]*)\)", parsed)

    # data are mapped to a dictionary
    for i in range(1, length(names_idx))
        if length(names_idx[i]) > 0
            # extract names using idx and split to separate the `=` sign.
            raw_names = [parsed[i][idx] for idx in names_idx[i]]
            # extract the naif ids
            naif_idxs = findall.(r"(?<=BODY)(.*?)(?=_)", raw_names)
            # extract the property name
            prop_idxs = findall.(r"(?<=\d_)(.*?)(?==)", raw_names)

            # trasform naifids to integers
            naif = parse.(Int64, [raw_names[j][ids[1]] for (j, ids) in enumerate(naif_idxs)])
            # trasform property names to symbols
            prop = Symbol.(lowercase.(strip.([raw_names[j][ids[1]] for (j, ids) in enumerate(prop_idxs)])))
            data = split.([replace(parsed[i][idx], "D" => "E") for idx in datas_idx[i]])

            for (name, body, value_) in zip(prop, naif, data)
                # parse a vector of floats
                value = parse.(Float64, value_[3:end-1])
                # create a temporary dictionary
                temp = Dict(body => Dict(name => length(value) > 1 ? value : value[1]))
                # merge with the global dictionary
                mergewith!(merge!, dict, temp)
            end
        end
    end 
    nothing
end