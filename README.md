# SMLMSim

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSMLM.github.io/SMLMSim.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSMLM.github.io/SMLMSim.jl/dev)
[![Build Status](https://github.com/JuliaSMLM/SMLMSim.jl/workflows/CI/badge.svg)](https://github.com/JuliaSMLM/SMLMSim.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaSMLM/SMLMSim.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaSMLM/SMLMSim.jl)

## Overview

SMLMSim is a Julia package for simulating Single Molecule Localization Microscopy (SMLM) data with realistic physical properties. It builds upon [SMLMData.jl](https://github.com/JuliaSMLM/SMLMData.jl), reexporting essential types and functions, and utilizes [MicroscopePSFs.jl](https://github.com/JuliaSMLM/MicroscopePSFs.jl) for realistic image generation.

The package provides tools for:

-   **Static SMLM Simulation:** Generating fixed spatial patterns (2D/3D) with realistic fluorophore photophysics (blinking) and localization uncertainty. Ideal for super-resolution studies.
-   **Diffusion & Interaction Simulation:** Modeling dynamic molecule behavior, including Brownian motion, dimerization, and dissociation using Smoluchowski dynamics. Suitable for single-particle tracking (SPT) studies.
-   **Microscope Image Generation:** Creating simulated camera images from emitter data using configurable Point Spread Functions (PSFs).

All simulations use physical units (microns, seconds) and produce data compatible with the broader [JuliaSMLM](https://github.com/JuliaSMLM) ecosystem.

## Installation

```julia
using Pkg
Pkg.add("SMLMSim")
```

## Quick Start

### Static SMLM Simulation

Simulate fixed patterns with blinking and localization noise.

```julia
using SMLMSim

# Define a camera and simulation parameters
camera = IdealCamera(128, 128, 0.1)  # 128×128 pixels, 100nm pixels
params = StaticSMLMParams(density=1.0, σ_psf=0.13) # Density 1/μm², PSF 130nm

# Run simulation for an 8-molecule ring pattern
smld_true, smld_model, smld_noisy = simulate(
    params; # Use semicolon to separate positional and keyword arguments
    pattern=Nmer2D(n=8, d=0.1), # 100nm diameter ring
    camera=camera
)

# smld_noisy contains realistic SMLM coordinates
println("Generated $(length(smld_noisy.emitters)) localizations.")
```
*Output:* `smld_true` (ground truth), `smld_model` (kinetics), `smld_noisy` (kinetics + noise).

### Diffusion & Interaction Simulation

Simulate molecules diffusing and interacting (e.g., dimerization).

```julia
using SMLMSim

# Set diffusion simulation parameters
params = DiffusionSMLMParams(
    density = 0.5,        # molecules per μm²
    box_size = 10.0,      # μm
    diff_monomer = 0.1,   # μm²/s
    k_off = 0.2,          # s⁻¹ dimer dissociation rate
    dt = 0.01,            # s simulation timestep
    t_max = 10.0          # s total simulation time
)

# Run diffusion simulation
smld = simulate(params) # Returns a BasicSMLD object with all emitters

println("Simulated diffusion for $(params.t_max) seconds.")
# 'smld' can be used for analysis or image generation
```

## Core Concepts

### Patterns

Define spatial arrangements (see `Pattern` types like `Nmer2D`, `Line3D`, `uniform2D`).

```julia
# Examples:
nmer = Nmer2D(n=8, d=0.1)  # 8 molecules in a 100nm diameter circle
line = Line3D(λ=5.0, endpoints=[(-1.0, 0.0, -0.5), (1.0, 0.0, 0.5)]) # 5 mols/μm
```

### Molecules & Photophysics

Model fluorophore behavior (e.g., `GenericFluor` with state transitions).

```julia
# Example: Two-state blinking model using positional constructor
fluor = GenericFluor(10000.0, [-10.0 10.0; 1e-2 -1e-2]) # γ=1e4, k_off=10, k_on=1e-2
```

### Localization Uncertainty

Realistic noise based on PSF width (`σ_psf`) and photon counts is added in static simulations.

### Image Generation

Create camera images from simulation results.

```julia
using MicroscopePSFs # Needed for PSF types

# Generate images from diffusion simulation output
psf = GaussianPSF(0.15) # 150nm PSF width
# Use smld_model to avoid double-counting localization errors
images = gen_images(smld, psf; 
    frame_integration=10, # 10 simulation time steps for each camera frame
    support=1.0 # PSF support range
    ) 

println("Generated $(size(images,3)) camera images.")
```

### sCMOS Camera with Realistic Noise

SMLMSim.jl 0.4+ supports realistic sCMOS camera noise modeling through SMLMData 0.4.

```julia
using SMLMSim
using MicroscopePSFs

# Create an sCMOS camera (128×128 pixels, 100nm pixels, 1.6 e⁻ read noise)
camera_scmos = SCMOSCamera(128, 128, 0.1, 1.6)

# Run static simulation with sCMOS camera
params = StaticSMLMParams(density=1.0, σ_psf=0.13)
smld_true, smld_model, smld_noisy = simulate(
    params,
    pattern=Nmer2D(n=8, d=0.1),
    camera=camera_scmos
)

# Generate images with full sCMOS noise model
# (quantum efficiency, Poisson, read noise, gain, offset)
psf = GaussianPSF(0.15)
images_scmos = gen_images(smld_noisy, psf, bg=10.0, camera_noise=true)

# For diffusion simulations
diff_params = DiffusionSMLMParams(density=0.5, box_size=10.0)
smld_diff = simulate(diff_params; camera=camera_scmos, override_count=10)
```

The sCMOS noise model applies:
1. **Quantum efficiency**: Photon → photoelectron conversion
2. **Poisson noise**: Shot noise on photoelectrons
3. **Read noise**: Gaussian noise per pixel
4. **Gain**: Electron → ADU conversion
5. **Offset**: Dark level addition

## Example Workflow: Static Simulation & Visualization

```julia
using SMLMSim
using CairoMakie # Requires installation: Pkg.add("CairoMakie")
using MicroscopePSFs

# --- Simulation Setup ---
camera = IdealCamera(128, 128, 0.1) # 128×128 pixels, 100nm pixels
params = StaticSMLMParams(density=1.0, σ_psf=0.13)
smld_true, smld_model, smld_noisy = simulate(
    params,
    pattern=Nmer2D(n=6, d=0.2), # Hexamer
    camera=camera
)

# --- Visualization ---
emitters = smld_noisy.emitters
x_coords = [e.x for e in emitters]
y_coords = [e.y for e in emitters]
photons = [e.photons for e in emitters]

fig = Figure(size=(600, 500))
ax = Axis(fig[1, 1],
    title="Simulated SMLM Localizations (Hexamer)",
    xlabel="x (μm)", ylabel="y (μm)",
    aspect=DataAspect(), yreversed=true
)
scatter!(ax, x_coords, y_coords, color=photons, colormap=:viridis, markersize=4, alpha=0.7)
Colorbar(fig[1, 2], colormap=:viridis, label="Photons")
display(fig)
# save("smlm_hexamer.png", fig)
```

## Further Information

For more detailed examples, API documentation, and explanations of the underlying models, please see the [Full Documentation](https://JuliaSMLM.github.io/SMLMSim.jl/dev).

## Contributors

-   [JuliaSMLM Team](https://github.com/JuliaSMLM)

## License

This project is licensed under the MIT License - see the LICENSE file for details.