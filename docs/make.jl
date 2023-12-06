using Documenter, FrameTransformations
using Pkg 

const CI = get(ENV, "CI", "false") == "true"

if CI 
    Pkg.add("Ephemerides")
    Pkg.add("ReferenceFrameRotations")
    Pkg.add("JSMDUtils")
end

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
            "04 - Loading EOP Data" => "Tutorials/t03_eop.md",
            "05 - Light Time Corrections" => "Tutorials/t04_lighttime.md",
            "06 - Multithreading" => "Tutorials/t05_multithread.md"
        ],
        
        "Use Case Examples" => [
            "CR3BP" => "Examples/e01_cr3bp.md",
            "High-Fidelity Earth-Moon Environment" => "Examples/e02_hifi.md"
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