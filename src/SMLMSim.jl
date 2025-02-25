module SMLMSim

using SMLMData
using Distributions
using LinearAlgebra

# Re-export critical types from SMLMData to make them available to users
export AbstractCamera, IdealCamera, AbstractEmitter, Emitter2D, Emitter3D, Emitter2DFit, Emitter3DFit, BasicSMLD

include("molecules.jl")
include("patterns.jl")
include("sim.jl")
include("interface.jl")

# Submodules
include("diffusion/InteractionDiffusion.jl")

# Import specific functions from InteractionDiffusion
using .InteractionDiffusion: DiffusingMolecule, DiffusingMoleculeSystem, 
                            SmoluchowskiParams, get_dimers, 
                            show_frame, visualize_sequence, visualize_simulation,
                            gen_image, gen_image_sequence

# Add this line to import the simulate method
using .InteractionDiffusion: simulate

# Export molecule types
export
    # Abstract molecule types
    Molecule,
    
    # Concrete molecule types
    GenericFluor

# Export simulation functions
export
    # Kinetic simulation
    intensitytrace,
    kineticmodel,
    noise,
    
    # CTMC type and functions
    CTMC,
    getstate,
    getnext

# Core types and functions for diffusion simulation
export
    # Core types
    DiffusingMolecule,
    DiffusingMoleculeSystem,
    SmoluchowskiParams,

    # Core simulation functions
    simulate, # Will dispatch to appropriate method based on argument types
    
    # Analysis functions
    get_dimers,
    gen_dimer_images,
    gen_image,         # Add this
    gen_image_sequence # Add this

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