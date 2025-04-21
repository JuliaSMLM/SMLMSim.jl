"""
    InteractionDiffusion

This module provides simulation tools for diffusion and interaction between particles.

# Overview
Simulates diffusion and interaction dynamics between particles in 2D/3D space.
Includes functionality for generating microscope images, analyzing dimers, and
visualizing particle dynamics.

# Components
- Abstract and concrete emitter types (AbstractDiffusingEmitter, DiffusingEmitter2D, DiffusingEmitter3D)
- Smoluchowski dynamics simulation
- Analysis tools for dimers and state transitions
- SMLD conversion utilities

# Examples
```julia
# Set up simulation parameters
params = DiffusionSMLMParams(
    density = 0.5,            # molecules per μm²
    box_size = 10.0,          # μm
    diff_monomer = 0.1,       # μm²/s
    diff_dimer = 0.05,        # μm²/s
    k_off = 0.2,              # s⁻¹
    r_react = 0.01,           # μm
    d_dimer = 0.05,           # μm
    dt = 0.01,                # s
    t_max = 10.0,             # s
    camera_framerate = 20.0,  # fps
    camera_exposure = 0.04    # s
)

# Run simulation - returns a single SMLD with all emitters
smld = simulate(params)

# Generate images for microscopy
psf = GaussianPSF(0.15)  # 150nm PSF width
images = gen_images(psf, smld)

# Analyze results
dimer_smld = get_dimers(smld)
frames, dimer_fractions = analyze_dimer_fraction(smld)
```
"""
module InteractionDiffusion

using SMLMData 
using Distributions
using MicroscopePSFs
using Statistics

using Printf

# Import the main simulate function to add our method
import ..simulate
import ..Core: AbstractSim
import ..Core: SMLMSimParams

include("types.jl")
include("smoluchowski.jl")
include("helpers.jl")
include("analysis.jl")

# Core types and functions for diffusion simulation
export
    # Core simulation parameters
    DiffusionSMLMParams,
    
    # Emitter types
    AbstractDiffusingEmitter,
    DiffusingEmitter2D,
    DiffusingEmitter3D,

    # Core simulation functions
    simulate,              # Returns a BasicSMLD with diffusing emitters
    
    # Analysis functions
    get_dimers,
    get_monomers,
    analyze_dimer_fraction,
    analyze_dimer_lifetime,
    track_state_changes,
    filter_by_state,
    
    # SMLD conversion utilities
    create_smld,
    get_frame

end