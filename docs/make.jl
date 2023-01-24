using Documenter, Basic

setup = quote
    using Basic
end
DocMeta.setdocmeta!(Basic, :DocTestSetup, setup; recursive = true)

makedocs(;
    authors="Astronaut Development Team",
    repo="https://gitlab.com/astronaut-tools/julia/core/Basic/blob/{commit}{path}#{line}",
    sitename="Basic",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://astronaut-tools.gitlab.io/julia/core/Basic",
        highlights = ["yaml"],
        ansicolor = true,
        assets=["assets/init.js"]
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => [
            "Epochs" => "Tutorials/t_01_epochs.md",
            "Timescales" => "Tutorials/t_02_extending_scales.md",
            "Ephemeris" => "Tutorials/t_03_ephem.md"
              ],
        # "Manual" => [
        #     "Overview" => "manual.md"
        # ],
        "Modules" => [
            "Tempo" => "Modules/time.md",
            "Ephemeris" => "Modules/ephem.md",
            "Orient" => "Modules/orient.md",
            "Frames" => "Modules/frames.md",
            "Utils" => "Modules/utils.md",
            "Graphs" => "Modules/graph.md"
        ]
    ],
    strict = !("strict=false" in ARGS),
    doctest = ("doctest=only" in ARGS) ? :only : true,
    # clean = false,
)
