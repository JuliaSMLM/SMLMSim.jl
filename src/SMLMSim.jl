"""
    SMLMSim

Main module for the SMLMSim.jl package.

# API Overview
For a comprehensive overview of the API, use the help mode on `api`:

    ?api

Or access the complete API documentation programmatically:

    docs = SMLMSim.api()

This package provides tools for simulating Single Molecule Localization Microscopy (SMLM) data.
It includes modules for:

- **Core:** Fundamental types (molecules, patterns) and photophysics simulation (CTMC, blinking).
- **StaticSMLM:** Simulating static emitters with blinking and localization noise.
- **InteractionDiffusion:** Simulating diffusing and interacting emitters (e.g., dimerization) using Smoluchowski dynamics.
- **CameraImages:** Generating simulated camera images from emitter data, including noise models.

The main `SMLMSim` module re-exports key types and functions from these submodules
to provide a unified user interface.

# Return Convention
All main functions return `(result, info)` tuples:
- `simulate()` returns `(smld, SimInfo)` - primary data plus metadata
- `gen_images()` returns `(images, ImageInfo)` - images plus metadata

# Usage
```julia
using SMLMSim

# Example: Static simulation - returns (smld_noisy, info)
params_static = StaticSMLMParams(density=1.0, σ_psf=0.13)
smld_noisy, info = simulate(params_static)
smld_true = info.smld_true    # Ground truth from info
smld_model = info.smld_model  # Kinetic model from info

# Example: Diffusion simulation - returns (smld, info)
params_diff = DiffusionSMLMParams(density=0.5, diff_monomer=0.1)
smld_diff, info = simulate(params_diff)

# Example: Generate images - returns (images, info)
psf = GaussianPSF(0.15)
images, img_info = gen_images(smld_noisy, psf)
```
"""
module SMLMSim

using SMLMData
using Distributions
using LinearAlgebra

# Re-export critical types from SMLMData to make them available to users
export AbstractCamera, IdealCamera, SCMOSCamera, AbstractEmitter, Emitter2D, Emitter3D, Emitter2DFit, Emitter3DFit, BasicSMLD

# Core module (includes molecules.jl and patterns.jl internally)
include("core/Core.jl")

include("interface.jl")

# Import specific functions from Core
using .Core: CTMC, get_state, get_next, intensity_trace, kinetic_model
using .Core: Molecule, GenericFluor, Pattern, Pattern2D, Pattern3D
using .Core: Nmer2D, Nmer3D, Line2D, Line3D, uniform2D, uniform3D, rotate!
using .Core: AbstractSim, SMLMSimParams # Add abstract types import
using .Core: get_track, get_num_tracks, get_tracks # Track utility functions
using .Core: SimInfo, ImageInfo # Info types for (result, info) return convention

# Include submodules after the Core imports are available
include("static/StaticSMLM.jl")
include("diffusion/InteractionDiffusion.jl")
include("camera_images/CameraImages.jl")

# Include the API overview functionality
include("api.jl")

# Import specific functions from InteractionDiffusion
using .InteractionDiffusion: DiffusionSMLMParams, get_dimers, 
                            get_monomers, analyze_dimer_fraction, analyze_dimer_lifetime,
                            DiffusingEmitter2D, DiffusingEmitter3D, extract_final_state

# Import from StaticSMLM
using .StaticSMLM: StaticSMLMParams, apply_noise

# Import from CameraImages
using .CameraImages: gen_images, gen_image, poisson_noise, poisson_noise!, scmos_noise, scmos_noise!

# Add this line to import the simulate methods
using .InteractionDiffusion: simulate
using .StaticSMLM: simulate

# Export simulation interfaces
export
    # Simulation interfaces
    AbstractSim,
    SMLMSimParams  # Add this export

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
    DiffusionSMLMParams,  # Changed here
    
    # New diffusing emitter types for imaging
    DiffusingEmitter2D,
    DiffusingEmitter3D,

    # Core simulation functions
    simulate, # Will dispatch to appropriate method based on argument types
    
    # Analysis functions
    get_dimers,
    get_monomers,
    analyze_dimer_fraction,
    analyze_dimer_lifetime,
    extract_final_state

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
    rotate!

# Export molecule types
export
    # Abstract molecule types
    Molecule,
    
    # Concrete molecule types
    GenericFluor,
    
    # Static SMLM types
    StaticSMLMParams,
    apply_noise

# Track utility functions
export
    get_track,
    get_num_tracks,
    get_tracks

# Visualization and imaging
export
    gen_images,
    gen_image,
    poisson_noise,
    poisson_noise!,
    scmos_noise,
    scmos_noise!,

# Info types for (result, info) return convention
    SimInfo,
    ImageInfo,

# API overview
    api

end