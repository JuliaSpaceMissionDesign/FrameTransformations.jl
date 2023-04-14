using Documenter, Basic

setup = quote
    using ReferenceFrameRotations
    using Basic
    using Basic.Tempo
    using Basic.Ephemeris 
    using Basic.Orient 
    using Basic.Frames 
    using Basic.Utils

    const DE440_KERNEL_PATH = "/home/andrea/Documents/Kernels/spk/de440.bsp"
end

DocMeta.setdocmeta!(Basic, :DocTestSetup, setup; recursive=true)

const CI = get(ENV, "CI", "false") == "true"

makedocs(;
    authors="Astronaut Development Team",
    repo="https://gitlab.com/astronaut-tools/julia/core/Basic/blob/{commit}{path}#{line}",
    sitename="Basic",
    modules=[Basic],
    format=Documenter.HTML(;
        prettyurls=CI,
        canonical="https://astronaut-tools.gitlab.io/julia/core/Basic",
        highlights=["yaml"],
        ansicolor=true,
        assets=["assets/init.js"],
    ),
    pages=[
        "Home" => "index.md",
        "Contributing" => "todos.md",
        "Tutorials" => [
            "Epochs" => "Tutorials/t_01_epochs.md",
            "Timescales" => "Tutorials/t_02_extending_scales.md",
            "Ephemeris" => "Tutorials/t_03_ephem.md",
            "Points Graphs" => "Tutorials/t_04_points.md",
            "Axes Graphs" => "Tutorials/t_05_axes.md",
            "Frames" => "Tutorials/t_06_frames.md",
        ],
        "Modules" => [
            "Tempo" => "Modules/time.md",
            "Ephemeris" => "Modules/ephem.md",
            "Orient" => "Modules/orient.md",
            "Frames" => "Modules/frames.md",
            "Utils" => "Modules/utils.md",
        ],
    ],
    # strict=!("strict=false" in ARGS),
    # doctest=("doctest=only" in ARGS) ? :only : true,
    clean = true,
)
