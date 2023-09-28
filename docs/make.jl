using Documenter, FrameTransformations
using Pkg 

const CREATE_TUTORIALS = true;

const CI = get(ENV, "CI", "false") == "true"

if CI

    Pkg.add("IJulia")
    Pkg.add("Conda") 
    using Conda, IJulia
    Conda.add("nbconvert")

    Pkg.add("ReferenceFrameRotations")
    using ReferenceFrameRotations

    Pkg.add("Ephemerides")
    using Ephemerides

    Pkg.add("CalcephEphemeris")
    using CalcephEphemeris

    const TUTORIAL_PATH = "docs/src/Tutorials"
    files = readdir(TUTORIAL_PATH)
    ipynb_files = filter(file -> endswith(file, ".ipynb"), files)

    # Convert
    for file in ipynb_files
        nbconvert = IJulia.find_jupyter_subcommand("nbconvert");
        append!(nbconvert.exec, ["--to", "markdown", "--execute", joinpath(TUTORIAL_PATH, file) ])
        run(nbconvert)
    end

end


makedocs(;
    authors="Julia Space Mission Design Development Team",
    sitename="FrameTransformations.jl",
    modules=[FrameTransformations],
    format=Documenter.HTML(; prettyurls=CI, highlights=["yaml"], ansicolor=true),
    pages=[
        "Home" => "index.md",
        "Modules" => [
            "Orient" => "Modules/orient.md",
            "Frames" => "Modules/frames.md",
            "Utils" => "Modules/utils.md",
        ],
        "Tutorials" => [
            "Points Graphs" => "Tutorials/t01_points.md",
            "Axes Graphs" => "Tutorials/t02_axes.md",
            "Frames Graphs" => "Tutorials/t03_frames.md"
        ]
    ],
    clean=true,
)

deploydocs(;
    repo="github.com/JuliaSpaceMissionDesign/FrameTransformations.jl", branch="gh-pages"
)