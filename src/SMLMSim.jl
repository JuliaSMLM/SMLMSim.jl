module SMLMSim

using SMLMData
using Distributions
using LinearAlgebra

# Re-export critical types from SMLMData to make them available to users
export AbstractCamera, IdealCamera, AbstractEmitter, Emitter2D, Emitter3D, Emitter2DFit, Emitter3DFit, BasicSMLD

# Core module (includes molecules.jl and patterns.jl internally)
include("core/Core.jl")

include("interface.jl")

# Import specific functions from Core
using .Core: CTMC, get_state, get_next, intensity_trace, kinetic_model
using .Core: Molecule, GenericFluor, Pattern, Pattern2D, Pattern3D
using .Core: Nmer2D, Nmer3D, Line2D, Line3D, uniform2D, uniform3D, rotate!

# Include submodules after the Core imports are available
include("static/StaticSMLM.jl")
include("diffusion/InteractionDiffusion.jl")
include("camera_images/CameraImages.jl")

# Import specific functions from InteractionDiffusion
using .InteractionDiffusion: DiffusingMolecule, DiffusingMoleculeSystem, 
                            SmoluchowskiParams, get_dimers, 
                            show_frame, visualize_sequence, visualize_simulation,
                            gen_image, gen_image_sequence, 
                            DiffusingEmitter2D, DiffusingEmitter3D

# Import from StaticSMLM
using .StaticSMLM: StaticSMLMParams, apply_noise

# Import from CameraImages
using .CameraImages: gen_images, gen_image

# Add this line to import the simulate methods
using .InteractionDiffusion: simulate
using .StaticSMLM: simulate

# Export simulation functions
export
    # Kinetic simulation
    intensity_trace, # Renamed from intensitytrace to match function name
    kinetic_model,   # Renamed from kineticmodel to match function name
    noise,
    
    # CTMC type and functions
    CTMC,
    get_state,
    get_next

# Core types and functions for diffusion simulation
export
    # Core types
    DiffusingMolecule,
    DiffusingMoleculeSystem,
    SmoluchowskiParams,
    
    # New diffusing emitter types for imaging
    DiffusingEmitter2D,
    DiffusingEmitter3D,

    # Core simulation functions
    simulate, # Will dispatch to appropriate method based on argument types
    
    # Analysis functions
    get_dimers,
    gen_dimer_images,
    gen_image,         
    gen_image_sequence 

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

# Export molecule types
export
    # Abstract molecule types
    Molecule,
    
    # Concrete molecule types
    GenericFluor,
    
    # Static SMLM types
    StaticSMLMParams,
    apply_noise

# Visualization and imaging
export
    show_frame,
    visualize_sequence,
    visualize_simulation,
    gen_images,
    gen_image

end