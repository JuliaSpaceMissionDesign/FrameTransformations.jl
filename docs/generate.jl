using Literate
using Dates

# TODO: Remove items from `SKIPFILE` as soon as they run on the latest stable 
ONLYSTATIC = []
EXAMPLE_DIRS = ["Tutorials", "Examples"]


function update_date(content)
    content = replace(content, "DATEOFTODAY" => Dates.DateTime(now()))
    return content
end

for edir in EXAMPLE_DIRS 

    gen_dir = joinpath(@__DIR__, "src", edir, "gen")
    example_dir = joinpath(@__DIR__, "src", edir)

    for example in filter!(x -> endswith(x, ".jl"), readdir(example_dir))
        input = abspath(joinpath(example_dir, example))
        script = Literate.script(input, gen_dir)
        code = strip(read(script, String))
        mdpost(str) = replace(str, "@__CODE__" => code)
        Literate.markdown(
            input, gen_dir, 
            preprocess = update_date,
            postprocess = mdpost,
            documenter = !(example in ONLYSTATIC))
        # Literate.notebook(input, gen_dir, execute = false)
    end

end

