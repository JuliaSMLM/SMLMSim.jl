using SMLMSim
using Documenter
using CairoMakie
using MicroscopePSFs
using Distributions

DocMeta.setdocmeta!(SMLMSim, :DocTestSetup, :(using SMLMSim, CairoMakie, MicroscopePSFs, Distributions); recursive=true)

makedocs(;
    modules=[SMLMSim],
    authors="klidke@unm.edu",
    sitename="SMLMSim.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaSMLM.github.io/SMLMSim.jl",
        assets=String[],
    ),
    doctest = false,
    pages=[
        "Home" => "index.md",
        "Core Components" => [
            "Patterns" => "core/patterns.md",
            "Photophysics" => "core/photophysics.md",
            "Localization Uncertainty" => "core/noise.md"
        ],
        "Static SMLM" => [
            "Overview" => "static/overview.md",
            "Examples" => "static/examples.md"
        ],
        "Diffusion-Interaction" => [
            "Overview" => "diffusion/overview.md",
            "Examples" => "diffusion/examples.md"
        ],
        "Microscope Images" => "images.md",
        "API Reference" => "api.md",
    ],
    warnonly=true,
)

deploydocs(;
    repo="github.com/JuliaSMLM/SMLMSim.jl",
    devbranch = "main"
)
