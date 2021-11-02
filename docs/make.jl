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
    ],
)

deploydocs(;
    repo="github.com/JuliaSMLM/SMLMSim.jl",
    devbranch = "main"
)
