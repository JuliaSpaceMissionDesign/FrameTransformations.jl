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
        assets=[]
    ),
    pages=[
        "Home" => "index.md",
        "Modules" => [
            "Bodies" => [
                "API" => "Bodies/api.md"
            ],
            "Tempo" => [
                "API" => "Tempo/api.md"
            ],
            "Orient" => [
                "API" => "Orient/api.md"
            ],
            "Rotate" => [
                "API" => "Rotate/api.md"
            ],
            "Utils" => [
                "API" => "Utils/api.md"
            ]
        ]
    ],
    strict = !("strict=false" in ARGS),
    doctest = ("doctest=only" in ARGS) ? :only : true,
    clean = false,
)