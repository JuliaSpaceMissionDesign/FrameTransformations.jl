using Documenter, FrameTransformations
using Pkg 

const CI = get(ENV, "CI", "false") == "true"

makedocs(;
    authors="Julia Space Mission Design Development Team",
    sitename="FrameTransformations.jl",
    modules=[FrameTransformations],
    format=Documenter.HTML(; prettyurls=CI, highlights=["yaml"], ansicolor=true),
    pages=[
        "Home" => "index.md",

        "Tutorials" => [
            "01 - Frame System" => "Tutorials/t00_frames.md",
            "02 - Axes" => "Tutorials/t01_axes.md",
            "03 - Points" => "Tutorials/t02_points.md",
            "04 - Use Case: CR3BP" => "Tutorials/t03_cr3bp.md",
            "05 - Use Case: High Fidelity" => "Tutorials/t04_hifi.md",
            "06 - Multithreading" => "Tutorials/t05_multithread.md"
        ],

        "Benchmarks" => "benchmarks.md",

        "Modules" => [
            "Frames" => [
                "Public API" => "Modules/frames_api.md",
                "Low-level API" => "Modules/frames_lapi.md",             
            ],

            "Orient" => [
                "Public API" => "Modules/orient_api.md"
                "Low-level API" => "Modules/orient_lapi.md"
            ],

            "Utils" => [
                "Public API" => "Modules/utils_api.md"
                "Low-level API" => "Modules/utils_lapi.md"
            ],

        ], 

        "Roadmap" => "roadmap.md"

    ],
    clean=true,
)

deploydocs(;
    repo="github.com/JuliaSpaceMissionDesign/FrameTransformations.jl", branch="gh-pages"
)