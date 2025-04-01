"""
    InteractionDiffusion

This module provides simulation tools for diffusion and interaction between particles.

# Overview
Simulates diffusion and interaction dynamics between particles in 2D/3D space.
Includes functionality for generating microscope images, analyzing dimers, and
visualizing particle dynamics.

# Components
- Core simulation types (DiffusingMolecule, DiffusingMoleculeSystem)
- Smoluchowski dynamics simulation
- Visualization tools
- Microscope image generation
- Dimer analysis and extraction

# Examples
```julia
# Set up simulation parameters
params = SmoluchowskiParams(
    density = 0.5,        # molecules per μm²
    box_size = 10.0,      # μm
    diff_monomer = 0.1,   # μm²/s
    diff_dimer = 0.05,    # μm²/s
    k_off = 0.2,          # s⁻¹
    r_react = 0.01,       # μm
    d_dimer = 0.05,       # μm
    dt = 0.01,            # s
    t_max = 10.0          # s
)

# Run simulation
systems = simulate(params)

# Visualize the simulation
visualize_sequence(systems, filename="diffusion.mp4")
```
"""
module InteractionDiffusion

using SMLMData 
using Distributions
using MicroscopePSFs

using Printf
using CairoMakie

# Import the main simulate function to add our method
import ..simulate
import ..SMLMSim: AbstractSim

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
    simulate,  # Export the function, but the implementation is a method
    simulate_and_image,

    # Analysis functions
    get_dimers,
    gen_dimer_images,
    analyze_dimer_fraction,

    # Visualization functions
    show_frame,
    visualize_sequence,
    visualize_simulation,
    
    # Microscope image generation
    gen_image,
    gen_image_sequence

end