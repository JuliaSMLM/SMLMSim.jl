```@meta
CurrentModule = SMLMSim
```

# SMLMSim.jl

*A Julia package for simulating Single Molecule Localization Microscopy (SMLM) data with realistic physical properties.*

## Overview

SMLMSim provides tools for generating single molecule localization microscopy data with physically realistic properties, including:

- Customizable spatial patterns of fluorophores in 2D and 3D
- Realistic fluorophore photophysics with stochastic kinetic models
- Accurate localization uncertainty based on photon counts
- Diffusion and interactions between molecules
- Microscope image generation with configurable PSFs

All simulations use physical units, with coordinates in microns and time in seconds, allowing for direct comparison with experimental data.

## Simulation Types

SMLMSim currently supports two main types of simulations:

### Static SMLM

Static simulations generate fixed patterns of molecules (such as protein complexes or structures) with realistic blinking behavior and localization uncertainty. This approach is ideal for SMLM super-resolution applications including:

- Simulating structured samples (oligomers, filaments, etc.)
- Testing localization algorithms
- Evaluating super-resolution reconstruction methods
- Benchmarking SMLM analysis software

### Diffusion-Interaction

Diffusion simulations model the dynamic behavior of molecules undergoing Brownian motion, suitable for single particle tracking applications, including:

- Free diffusion with configurable coefficients
- Formation of molecular complexes (dimerization)
- Dissociation of complexes
- Combined translational and rotational diffusion

This approach is ideal for studying dynamic biological processes and single-particle tracking applications.

## Installation

```julia
using Pkg
Pkg.add("SMLMSim")
```

## Quick Start

### Static SMLM Simulation

```julia
using SMLMSim

# Define a camera with 100nm pixel size
camera = IdealCamera(128, 128, 0.1)

# Run a basic static simulation
smld_true, smld_model, smld_noisy = simulate(
    density=1.0,                # 1 pattern per μm²
    σ_psf=0.13,           # 130nm PSF width
    pattern=Nmer2D(n=8, d=0.1),  # 8-molecule circular pattern (100nm diameter)
    camera=camera
)
```

### Diffusion-Interaction Simulation

```julia
# Set diffusion simulation parameters
params = DiffusionSMLMParams(
    density = 0.5,        # molecules per μm²
    box_size = 10.0,      # μm
    diff_monomer = 0.1,   # μm²/s
    diff_dimer = 0.05,    # μm²/s
    k_off = 0.2,          # s⁻¹
    dt = 0.01,            # s
    t_max = 10.0          # s
)

# Run diffusion simulation
smld_diffusion = simulate(params)
```

## Package Structure

SMLMSim is built upon [SMLMData.jl](https://github.com/JuliaSMLM/SMLMData.jl), reexporting essential types and functions so you typically don't need to import SMLMData directly.

The main components of the package are:

- **Patterns**: Spatial arrangements of molecules (Nmer2D, Line2D, etc.)
- **Molecules**: Photophysical models (e.g., GenericFluor)
- **Simulation**: Kinetic models and noise generation
- **Diffusion**: Smoluchowski dynamics for molecular interactions



## License

This project is licensed under the MIT License - see the LICENSE file for details.