using Documenter, FrameTransformations

const CI = get(ENV, "CI", "false") == "true"

makedocs(;
    authors="Andrea Pasquale, Michele Ceresoli and contributors",
    sitename="FrameTransformations.jl",
    modules=[FrameTransformations],
    format=Documenter.HTML(;
        prettyurls=CI,
        highlights=["yaml"],
        ansicolor=true,
    ),
    pages=[
        "Home" => "index.md",
        "Modules" => [
            "Orient" => "Modules/orient.md",
            "Frames" => "Modules/frames.md",
            "Utils" => "Modules/utils.md",
        ],
    ],
    clean=true,
)
