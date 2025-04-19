```@meta
CurrentModule = SMLMSim
```

# Diffusion-Interaction Simulation

## Overview

The Diffusion-Interaction module simulates dynamic molecular processes including diffusion and interactions between particles in a controlled environment. This allows you to model realistic molecular behaviors such as:

- Free diffusion of monomers
- Formation of molecular complexes (dimers)
- Dissociation of complexes
- Combined translational and rotational diffusion

The simulation operates within a defined box with customizable physical parameters, making it suitable for studying a wide range of biological phenomena at the single-molecule level.

## Simulation Model

The diffusion simulation is based on the Smoluchowski dynamics model with the following components:

1. **Particles**: Represented as point particles (monomers) or rigid structures (dimers)
2. **Diffusion**: Isotropic Brownian motion with specified diffusion coefficients
3. **Reactions**: 
   - Association: Two monomers within reaction radius form a dimer
   - Dissociation: Dimers break with rate k_off
4. **Boundaries**: Periodic or reflecting boundary conditions

At each time step, the simulation:
- Updates molecular states (dimerization/dissociation)
- Updates positions with appropriate diffusion models
- Handles boundary conditions

### Physical Units

All simulation parameters use consistent physical units:
- Spatial dimensions: microns (μm)
- Time: seconds (s)
- Diffusion coefficients: μm²/s
- Rate constants: s⁻¹

## Getting Started

### Running a Basic Simulation

The main interface for running diffusion simulations is the `simulate` function with `DiffusionSMLMParams`:

```julia
using SMLMSim

# Set simulation parameters
params = DiffusionSMLMParams(
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

# Run the simulation
smld = simulate(params)
```

The `smld` output is a `BasicSMLD` structure containing all emitters across all time points, with each emitter having frame information corresponding to the camera settings.

## Simulation Parameters

The `DiffusionSMLMParams` structure allows you to customize various aspects of the simulation:

```julia
# More complex simulation
params = DiffusionSMLMParams(
    density = 2.0,           # Higher density
    box_size = 20.0,         # Larger area
    diff_monomer = 0.2,      # Faster monomer diffusion
    diff_dimer = 0.08,       # Faster dimer diffusion
    diff_dimer_rot = 1.0,    # Faster rotational diffusion
    k_off = 0.05,            # Slower dissociation (more stable dimers)
    r_react = 0.015,         # Larger reaction radius
    d_dimer = 0.08,          # Larger dimer separation
    dt = 0.005,              # Smaller time step (higher precision)
    t_max = 30.0,            # Longer simulation
    ndims = 2,               # 2D simulation (default)
    boundary = "periodic",   # Periodic boundaries (default)
    camera_framerate = 20.0, # Camera frames per second
    camera_exposure = 0.04   # Camera exposure time
)
```



## Microscope Image Generation

The diffusion simulation can be converted into realistic microscope images using point spread function models.

### Creating Images from Simulation Data

```julia
# Set up camera and PSF
pixelsize = 0.1  # 100nm pixels
pixels = Int64(round(params.box_size/pixelsize))
camera = IdealCamera(1:pixels, 1:pixels, pixelsize)

# Set up PSF (Gaussian with 150nm width)
using MicroscopePSFs
psf = MicroscopePSFs.Gaussian2D(0.15)  # 150nm PSF width

# Generate images
image_stack = gen_images(smld, psf;
    photons=1000.0,
    bg=5.0,
    poisson_noise=true
)
```

## Analyzing Results

### Extracting Dimers

To focus on the behavior of dimers, you can extract only the molecules in dimer state:

```julia
# Extract dimers from the full simulation
dimer_smld = get_dimers(smld)
```

### Analyzing Dimer Formation

To quantify the formation of dimers over time:

```julia
# Calculate fraction of molecules in dimer state
frames, dimer_fractions = analyze_dimer_fraction(smld)

# Plot dimer formation over time
using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1],
    xlabel="Frame",
    ylabel="Fraction of molecules in dimers",
    title="Dimer formation dynamics"
)
lines!(ax, frames, dimer_fractions)
fig
```

## Emitter Types

The diffusion module introduces specialized emitter types:

- `DiffusingEmitter2D`: 2D emitter with state information (monomer/dimer)
- `DiffusingEmitter3D`: 3D emitter with state information (monomer/dimer)

These types include additional properties:
- `timestamp`: Actual simulation time
- `state`: Molecular state (`:monomer` or `:dimer`)
- `partner_id`: ID of linked molecule (for dimers)

