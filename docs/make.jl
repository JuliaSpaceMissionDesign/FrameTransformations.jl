using Documenter, FrameTransformations
using Pkg

const CI = get(ENV, "CI", "false") == "true"

if CI
    Pkg.add("Ephemerides")
    Pkg.add("StaticArrays")
    Pkg.add("ReferenceFrameRotations")
    Pkg.add("JSMDUtils")
    Pkg.add("JSMDInterfaces")
    Pkg.add("Literate")
    Pkg.add("Dates")
    Pkg.add("Tempo")

    # Examples 
    Pkg.add("DiffEqBase")
    Pkg.add("OrdinaryDiffEq")
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
            "04 - Points" => "Tutorials/gen/t03_points.md"
        ],
        "Examples" => [
            "01 - Custom Orbit" => "Examples/gen/e00_ode.md"
        ],
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
