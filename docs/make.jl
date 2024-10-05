using Documenter, FrameTransformations
using Pkg

const CI = get(ENV, "CI", "false") == "true"

if CI
    Pkg.add("Ephemerides")
    Pkg.add("ReferenceFrameRotations")
    Pkg.add("JSMDUtils")
    Pkg.add("JSMDInterfaces")
    Pkg.add("Literate")
    Pkg.add("Dates")
    Pkg.add("Tempo")
end

include("generate.jl")

makedocs(;
    authors="JSMD Development Team",
    sitename="FrameTransformations.jl",
    modules=[FrameTransformations],
    format=Documenter.HTML(; prettyurls=CI, highlights=["yaml"], ansicolor=true),
    pages=[
        "Home" => "index.md",
        "Tutorials" => [
            "01 - Frame System" => "Tutorials/gen/t00_frames.md",
            "02 - Rotation" => "Tutorials/gen/t01_rotation.md",
            "03 - Axes" => "Tutorials/gen/t02_axes.md",
            "04 - Points" => "Tutorials/gen/t03_points.md",
            # "05 - Light Time Corrections" => "Tutorials/gen/t04_lighttime.md",
            # "06 - Multi-threading" => "Tutorials/gen/t05_multithread.md"
        ],
        # "Use Cases" => [
        #     "CR3BP" => "Examples/gen/e01_cr3bp.md",
        #     "High Fidelity" => "Examples/gen/e02_hifi.md",
        #     "Custom Orbit" => "Examples/gen/e03_customorb.md"
        # ],
        # "Benchmarks" => "benchmarks.md",
        "API" => [
            "Public API" => [
                "Axes" => "API/axes_api.md",
                "Points" => "API/point_api.md",
                "Directions" => "API/dir_api.md",
                "Frames" => "API/frames_api.md"
            ],
        ],
    ],
    clean=true,
    checkdocs=:none
)

if CI
    deploydocs(;
        repo="github.com/JuliaSpaceMissionDesign/FrameTransformations.jl", branch="gh-pages"
    )
end
