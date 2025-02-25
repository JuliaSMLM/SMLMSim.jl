using SMLMSim
using Documenter

DocMeta.setdocmeta!(SMLMSim, :DocTestSetup, :(using SMLMSim); recursive=true)

makedocs(;
    modules=[SMLMSim],
    authors="klidke@unm.edu",
    repo="https://github.com/JuliaSMLM/SMLMSim.jl/blob/{commit}{path}#{line}",
    sitename="SMLMSim.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaSMLM.github.io/SMLMSim.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Core Concepts" => "core_concepts.md",
        "User Guide" => [
            "Basic Simulation" => "user_guide/basic_simulation.md",
            "Patterns" => "user_guide/patterns.md",
            "Photophysics" => "user_guide/photophysics.md",
            "Localization Uncertainty" => "user_guide/uncertainty.md",
        ],
        "Interaction-Diffusion" => "diffusion.md",
        "Examples" => "examples.md"
        # "Performance Tips" => "performance_tips.md",
        # "API Reference" => "api.md",
    ],
    warnonly=true,
)

deploydocs(;
    repo="github.com/JuliaSMLM/SMLMSim.jl",
    devbranch = "main"
)