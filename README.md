# SMLMSim

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSMLM.github.io/SMLMSim.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSMLM.github.io/SMLMSim.jl/dev)
[![Build Status](https://github.com/JuliaSMLM/SMLMSim.jl/workflows/CI/badge.svg)](https://github.com/JuliaSMLM/SMLMSim.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaSMLM/SMLMSim.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaSMLM/SMLMSim.jl)

## Overview

SMLMSim is a Julia package for simulating Single Molecule Localization Microscopy (SMLM) data with realistic physical properties. The package builds upon SMLMData, reexporting essential types and functions so you typically don't need to import SMLMData directly.

The package provides tools for:

- Generating spatial patterns of fluorophores in 2D and 3D
- Simulating fluorophore photophysics with stochastic kinetic models
- Adding realistic localization uncertainty based on photon counts
- Simulating diffusion and interactions between molecules
- Generating microscope images with configurable PSFs

All simulations use physical units, with coordinates in microns and time in seconds. The resulting data is organized into `SMLMData.SMLD` structures compatible with the broader JuliaSMLM ecosystem.

## Installation

```julia
using Pkg
Pkg.add("SMLMSim")
```

## Basic Usage

The high-level interface for simulating SMLM super-resolution coordinate data is the `simulate()` function (or the alias `sim()`).

```julia
using SMLMSim

# Basic simulation with default parameters
camera = IdealCamera(1:128, 1:128, 0.1)  # 128×128 pixels, 100nm pixels
smld_true, smld_model, smld_noisy = simulate(
    camera=camera
)
```

This basic example creates a 2D simulation using default parameters:
- 8-molecule circular patterns (Nmer2D with n=8, d=0.1μm)
- 1 pattern per square micron (ρ=1.0)
- PSF width of 130nm (σ_psf=0.13μm)
- Two-state fluorophore kinetics with realistic blinking behavior
- 1000 frames at 50 frames per second
- Minimum photon threshold of 50 for detection

The function returns three SMLD objects with SMLM coordinates:
- `smld_true`: Ground truth emitter positions (spatial coordinates only)
- `smld_model`: Positions with simulated blinking kinetics (subset of true positions appearing in different frames)
- `smld_noisy`: Positions with both blinking and localization uncertainty (realistic SMLM data with position errors)

For more control, you can customize the parameters:

```julia
# More customized simulation
smld_true, smld_model, smld_noisy = simulate(;
    ρ=1.0,                # emitters per μm²
    σ_psf=0.13,           # PSF width in μm (130nm)
    minphotons=50,        # minimum photons for detection
    ndatasets=10,         # number of independent datasets
    nframes=1000,         # frames per dataset
    framerate=50.0,       # frames per second
    pattern=Nmer2D(n=6, d=0.2),  # hexamer with 200nm diameter
    molecule=GenericFluor(; q=[0 50; 1e-2 0]),  # rates in 1/s
    camera=IdealCamera(; ypixels=256, xpixels=128, pixelsize=0.1)  # pixelsize in μm
)
```

This customized example:
- Creates hexagonal patterns (6 molecules in a 200nm circle)
- Uses 10 independent datasets of 1000 frames each
- Simulates fluorophores with specific on/off transition rates
- Uses a rectangular camera field of view (128×256 pixels)

In both cases, the output SMLD objects contain emitter information (x, y coordinates, photon counts, frame numbers, etc.) that can be used for further analysis or visualization.

## Pattern Types

SMLMSim includes several built-in pattern types for positioning fluorophores:

### 2D Patterns

```julia
# N molecules arranged in a circle
nmer = Nmer2D(n=8, d=0.1)  # 8 molecules in a 100nm diameter circle

# Linear pattern with random positions
line = Line2D(λ=5.0, endpoints=[(-2.0, 0.0), (2.0, 0.0)])  # 5 molecules per μm along line
```

### 3D Patterns

```julia
# N molecules arranged in a circle at z=0
nmer3d = Nmer3D(n=8, d=0.1)  # 8 molecules in a 100nm diameter circle

# 3D line with random positions
line3d = Line3D(λ=5.0, endpoints=[(-1.0, 0.0, -0.5), (1.0, 0.0, 0.5)])
```

## Molecule Models

SMLMSim supports different fluorophore photophysical models:

```julia
# Generic fluorophore with two-state kinetics
fluor = GenericFluor(
    γ=10000.0,           # photon emission rate in Hz
    q=[0 10; 1e-2 0]     # transition rate matrix: state 1 ↔ state 2
)
```

## Diffusion and Interaction Simulation

The package includes tools for simulating diffusion and interactions between molecules:

```julia
# Set up parameters for Smoluchowski diffusion simulation
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
visualize_sequence(systems, filename="diffusion.mp4", framerate=round(Int64,1/params.dt))

# Generate microscope images
psf = Gaussian2D(0.15)  # 150nm PSF width
images = gen_image_sequence(
    psf, 
    systems,
    frame_integration=10
)

# Extract only dimers
dimer_systems = get_dimers(systems)
dimer_images = gen_image_sequence(
    psf, 
    dimer_systems, 
    frame_integration=10
)
```

## Example Workflows

### 2D Simulation with Visualization

```julia
using SMLMSim
using CairoMakie

# Create camera with physical pixel size
camera = IdealCamera(1:128, 1:256, 0.1)  # 128×256 pixels, 100nm pixels

# Simulation parameters in physical units
smld_true, smld_model, smld_noisy = simulate(;
    ρ=1.0,                # emitters per μm²
    σ_psf=0.13,           # PSF width in μm
    pattern=Nmer2D(n=6, d=0.2),  # hexamer with 200nm diameter
    camera=camera
)

# Extract coordinates from emitters
x_noisy = [e.x for e in smld_noisy.emitters]
y_noisy = [e.y for e in smld_noisy.emitters]
photons = [e.photons for e in smld_noisy.emitters]

# Create figure and plot results
fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1], 
    title="Simulated SMLM Localizations",
    xlabel="x (μm)",
    ylabel="y (μm)",
    aspect=DataAspect(),
    yreversed=true  # This makes (0,0) at top-left
)

# Scatter plot with photon counts as color
scatter!(ax, x_noisy, y_noisy, 
    color=photons,
    colormap=:viridis,
    markersize=4,
    alpha=0.6
)

Colorbar(fig[1, 2], colormap=:viridis, label="Photons")

# Show or save the figure
display(fig)
# save("smlm_simulation.png", fig)
```

### 3D Simulation

```julia
using SMLMSim

# Create camera with physical pixel size
camera = IdealCamera(1:128, 1:256, 0.1)  # 128×256 pixels, 100nm pixels

# Simulation parameters in physical units
smld_true, smld_model, smld_noisy = simulate(;
    ρ=0.5,                # emitters per μm²
    pattern=Nmer3D(n=8, d=0.3),  # 3D pattern
    camera=camera,
    zrange=[-2.0, 2.0]    # 4μm axial range
)
```

## Contributors

- [JuliaSMLM Team](https://github.com/JuliaSMLM)

## License

This project is licensed under the MIT License - see the LICENSE file for details. The MIT License is a permissive license that allows for reuse with few restrictions. It permits use, modification, distribution, and private use while preserving copyright and license notices.