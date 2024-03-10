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

    Literate.markdown(joinpath("Tutorials", "t00_frames.jl"), joinpath("Tutorials"))
    Literate.markdown(joinpath("Tutorials", "t01_axes.jl"), joinpath("Tutorials"))
    Literate.markdown(joinpath("Tutorials", "t02_points.jl"), joinpath("Tutorials"))

end

# include("generate.jl")

makedocs(;
    authors="Julia Space Mission Design Development Team",
    sitename="FrameTransformations.jl",
    modules=[FrameTransformations],
    format=Documenter.HTML(; prettyurls=CI, highlights=["yaml"], ansicolor=true),
    pages=[

        "Home" => "index.md",

        "Tutorials" => [
            "01 - Frame System" => joinpath("Tutorials", "t00_frames.md"),
            "02 - Axes" => joinpath("Tutorials", "t01_axes.md"),
            "03 - Points" => joinpath("Tutorials", "t02_points.md"),
            # "04 - Loading EOP Data" => "Tutorials/gen/t03_eop.md",
            # "05 - Light Time Corrections" => "Tutorials/gen/t04_lighttime.md",
            # "06 - Multithreading" => "Tutorials/gen/t05_multithread.md"
        ],

        # "Use Cases" => [
        #     "CR3BP" => "Examples/gen/e01_cr3bp.md",
        #     "High-Fidelity Earth-Moon Environment" => "Examples/gen/e02_hifi.md",
        #     "Custom Orbit Representation" => "Examples/gen/e03_customorb.md"
        # ],

        # "Benchmarks" => "benchmarks.md",

        "API" => [
            "Public API" => "API/frames_api.md",
            "Low-level API" => "API/frames_lapi.md",             
        ],


        # "Roadmap" => "roadmap.md"

    ],
    clean=true,
)

deploydocs(;
    repo="github.com/JuliaSpaceMissionDesign/FrameTransformations.jl", branch="gh-pages"
)