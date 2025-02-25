```@meta
CurrentModule = SMLMSim
```

# SMLMSim.jl

*A Julia package for simulating Single Molecule Localization Microscopy (SMLM) data with realistic physical properties.*

## Overview

SMLMSim provides tools for generating super-resolution microscopy data with physically realistic properties, including:

- Customizable spatial patterns of fluorophores in 2D and 3D
- Realistic fluorophore photophysics with stochastic kinetic models
- Accurate localization uncertainty based on photon counts
- Diffusion and interactions between molecules
- Microscope image generation with configurable PSFs

All simulations use physical units, with coordinates in microns and time in seconds, allowing for direct comparison with experimental data.

## Installation

```julia
using Pkg
Pkg.add("SMLMSim")
```

## Quick Start

The high-level interface for simulating SMLM coordinate data is the `simulate()` function:

```julia
using SMLMSim

# Define a camera with 100nm pixel size
camera = IdealCamera(1:128, 1:128, 0.1)

# Run a basic simulation
smld_true, smld_model, smld_noisy = simulate(
    ρ=1.0,                # 1 pattern per μm²
    σ_psf=0.13,           # 130nm PSF width
    pattern=Nmer2D(n=8, d=0.1),  # 8-molecule circular pattern (100nm diameter)
    camera=camera
)
```

This generates:
- `smld_true`: Ground truth emitter positions
- `smld_model`: Positions with realistic blinking behavior
- `smld_noisy`: Positions with both blinking and localization uncertainty

## Package Structure

SMLMSim is built upon [SMLMData.jl](https://github.com/JuliaSMLM/SMLMData.jl), reexporting essential types and functions so you typically don't need to import SMLMData directly.

The main components of the package are:

- **Patterns**: Spatial arrangements of molecules (Nmer2D, Line2D, etc.)
- **Molecules**: Photophysical models (e.g., GenericFluor)
- **Simulation**: Kinetic models and noise generation
- **Interaction-Diffusion**: Smoluchowski dynamics for molecular interactions

## Citation

If you use this package in your research, please cite:

## Citing SMLMSim.jl

If you use SMLMSim.jl in your research, please cite this package as:

> Lidke, K. et al. (2025). SMLMSim.jl [Computer software]. https://github.com/JuliaSMLM/SMLMSim.jl

## Contributing

Contributions are welcome! Please see our [contributing guidelines](https://github.com/JuliaSMLM/SMLMSim.jl/blob/main/CONTRIBUTING.md) for more information.

```@contents
Pages = ["core_concepts.md", "diffusion.md", "examples.md", "performance_tips.md", "api.md"]
Depth = 2
```