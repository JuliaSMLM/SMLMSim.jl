"""
This module provides simulation tools for diffusion
    and interaction between particles.
"""
module InteractionDiffusion

using Distributions
using CairoMakie
using MicroscopePSFs
using SMLMSim

include("types.jl")
include("smoluchowski.jl")
include("visualize.jl")
include("microscope.jl")
include("dimer.jl")

# Core types and functions for diffusion simulation
export
    # Core types
    DiffusingMolecule,
    DiffusingMoleculeSystem,
    SmoluchowskiParams,

    # Core simulation functions
    simulate,
    simulate_and_image,

    # Analysis functions
    get_dimers,
    gen_dimer_images

end
