module SMLMSim

using SMLMData
using Distributions
using LinearAlgebra


include("molecules.jl")
include("patterns.jl")
include("sim.jl")
include("interface.jl")

# Submodules
include("diffusion/InteractionDiffusion.jl")

using SMLMSim.InteractionDiffusion


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

# Pattern simulation types and functions
export
    # Abstract pattern types
    Pattern,
    Pattern2D,
    Pattern3D,

    # Concrete patterns
    Nmer2D,
    Nmer3D,
    Line2D,
    Line3D,

    # Pattern generation
    uniform2D,
    uniform3D,
    rotate!,

    # Full simulation pipeline
    sim

# Visualization
export
    show_frame,
    visualize_sequence,
    visualize_simulation


end

