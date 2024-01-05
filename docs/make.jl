using Documenter, FrameTransformations
using Pkg 

const CI = get(ENV, "CI", "false") == "true"

if CI 
    Pkg.add("Ephemerides")
    Pkg.add("ReferenceFrameRotations")
    Pkg.add("JSMDUtils")
    Pkg.add("Literate")
    Pkg.add("Dates")
end

include("generate.jl")

makedocs(;
    authors="Julia Space Mission Design Development Team",
    sitename="FrameTransformations.jl",
    modules=[FrameTransformations],
    format=Documenter.HTML(; prettyurls=CI, highlights=["yaml"], ansicolor=true),
    pages=[
        "Home" => "index.md",

        "Tutorials" => [
            "01 - Frame System" => "Tutorials/gen/t00_frames.md",
            "02 - Axes" => "Tutorials/gen/t01_axes.md",
            "03 - Points" => "Tutorials/gen/t02_points.md",
            "04 - Loading EOP Data" => "Tutorials/gen/t03_eop.md",
            "05 - Light Time Corrections" => "Tutorials/gen/t04_lighttime.md",
            "06 - Multithreading" => "Tutorials/gen/t05_multithread.md"
        ],

        "Use Cases" => [
            "CR3BP" => "Examples/gen/e01_cr3bp.md",
            "High-Fidelity Earth-Moon Environment" => "Examples/gen/e02_hifi.md",
            "Custom Orbit Representation" => "Examples/gen/e03_customorb.md"
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

        ], 

        "Roadmap" => "roadmap.md"

    ],
    clean=true,
)

deploydocs(;
    repo="github.com/JuliaSpaceMissionDesign/FrameTransformations.jl", branch="gh-pages"
)