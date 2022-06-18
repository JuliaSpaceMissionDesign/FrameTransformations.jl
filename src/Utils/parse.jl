export parse_pck

"""
    parse_pck(filename:: String)
Open a NAIF ASCII PCK file and parse its data in a dictionary 
with Symbol => Vector{Float64} structure.
"""
function parse_pck(filename::String)
    # read and strip lines (remove tabs and spaces)
    lines = strip.(readlines(filename))
    # extract lines which are within `\begindata` and `\begintext`
    parsed = split(join(lines, " "), r"(?<=\\begintext).*?(?=\\begindata\s*BODY*)") 
    
    # extract lines which actually have data using the `BODY**** =` pattern
    # this is useful to filter the header and eventual trailings
    # this vector contains a list of `BODY******* =` elements which will be splitted afterwards
    names_idx = findall.(r"(BODY\w{1,}\D{1,}=)", parsed)
    # row data are extracted as between square brackets, the `=` before the brackets is added
    raw_datas_idx = findall.(r"=\D{1,}\(([^()]*)\)", parsed)

    # data are mapped to a dictionary
    mapped = Dict{Symbol, Vector{Float64}}()
    for i in range(1, length(names_idx))
        if length(names_idx[i]) > 0
            # extract names using idx and split to separate the `=` sign.
            sliced_names = split.([parsed[i][idx] for idx in names_idx[i]])
            # extract data and replace unsopported characters
            # split to separate the initial `= (` and the final `)`.
            sliced_datas = split.([replace(parsed[i][idx], "D" => "E") for idx in raw_datas_idx[i]])
            for j in range(1, length(sliced_names))
                # `parse` is called to transfrom strings into Float64
                push!(mapped, Symbol(sliced_names[j][1]) => parse.(Float64, sliced_datas[j][3:end-1]))
            end
        end
    end
    return mapped
end