```@meta
CurrentModule = SMLMSim.InteractionDiffusion
DocTestSetup = quote
    using SMLMSim
end
```

# Interaction-Diffusion

## Overview

The Interaction-Diffusion module simulates dynamic molecular processes including diffusion and interactions between particles in a controlled environment. This allows you to model realistic molecular behaviors such as:

- Free diffusion of monomers
- Formation of molecular complexes (dimers)
- Dissociation of complexes
- Combined translational and rotational diffusion

The simulation operates within a defined box with customizable physical parameters, making it suitable for studying a wide range of biological phenomena at the single-molecule level.

## Core Concepts

### Simulation Model

The diffusion simulation is based on the Smoluchowski dynamics model with the following components:

1. **Particles**: Represented as point particles (monomers) or rigid structures (dimers)
2. **Diffusion**: Isotropic Brownian motion with specified diffusion coefficients
3. **Reactions**: 
   - Association: Two monomers within reaction radius form a dimer
   - Dissociation: Dimers break with rate k_off
4. **Boundaries**: Periodic or reflecting boundary conditions

At each time step, the simulation:
- Updates the molecular states (dimerization/dissociation)
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

The main interface for running diffusion simulations is the `simulate` function with `SmoluchowskiParams`:

```julia
using SMLMSim

# Set simulation parameters
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

# Run the simulation
systems = simulate(params)
```

The `systems` output is a vector of `DiffusingMoleculeSystem` objects, each representing the state of the system at a specific time point. The time between frames is determined by the `dt` parameter.

### Customizing Simulation Parameters

The `SmoluchowskiParams` structure allows you to customize various aspects of the simulation:

```julia
# More complex simulation
params = SmoluchowskiParams(
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
    boundary = "periodic"    # Periodic boundaries (default)
)
```

### Parameter Guidelines

For realistic simulations, consider these guidelines:

| Parameter | Typical Range | Notes |
|-----------|---------------|-------|
| `density` | 0.1-10 μm⁻² | Depends on biological system |
| `diff_monomer` | 0.01-10 μm²/s | 0.1 μm²/s ≈ small protein in cytoplasm |
| `diff_dimer` | 0.5-0.8× monomer | Scales roughly with size |
| `k_off` | 0.001-10 s⁻¹ | Impacts complex stability |
| `r_react` | 0.001-0.02 μm | Reaction distance |
| `d_dimer` | 0.01-0.1 μm | Physical size of complex |

## Visualization

### Displaying a Single Frame

To visualize the state of the system at a specific time point:

```julia
# Show a specific frame
framenum = 50  # Time point index
show_frame(systems[framenum])
```

This produces a scatter plot with:
- Blue dots: Monomers
- Red dots: Dimers
- Position: Current location in the simulation box

For saving the visualization to a file:

```julia
# Save frame to file
show_frame(systems[framenum], "diffusion_frame50.png")
```

### Creating Animations

To create a movie showing the entire simulation:

```julia
# Generate an MP4 animation
visualize_sequence(systems, filename="diffusion_simulation.mp4")
```

This function renders each time point and compiles them into a movie, with optional parameters:
- `framerate`: Controls playback speed
- `show_dimers`: Whether to color-code dimers differently

## Microscope Image Generation

The diffusion simulation can be converted into realistic microscope images using point spread function models.

### Creating Images from a Single Frame

```julia
# Set up camera and PSF
pixelsize = 0.1  # 100nm pixels
pixels = Int64(round(params.box_size/pixelsize))
camera = IdealCamera(128, 128, pixelsize)

# Set up PSF (Gaussian with 150nm width)
using MicroscopePSFs
psf = MicroscopePSFs.Gaussian2D(0.15)  # 150nm PSF width

# Generate image for a specific frame
frame_image = gen_image(psf, systems[50], 1;
    photons=1000.0,  # photons per molecule
    bg=5.0,          # background photons per pixel
    poisson_noise=true
)
```

### Generating Image Sequences

For a complete movie with frame integration:

```julia
# Generate full image sequence
image_stack = gen_image_sequence(psf, systems;
    photons=1000.0,
    bg=5.0,
    frame_integration=10,  # integrate 10 simulation frames per output frame
    poisson_noise=true
)
```

The `frame_integration` parameter allows longer camera exposure times relative to the simulation time step, which is often realistic for experimental scenarios.

### One-Step Simulation and Imaging

For convenience, you can simulate and generate images in one step:

```julia
# Simulation and imaging in one step
images, systems = simulate_and_image(params, psf;
    photons=2000.0,
    bg=10.0,
    frame_integration=5
)
```

## Analyzing Dimers

### Extracting Dimers

To focus on the behavior of dimers, you can extract only the molecules in dimer state:

```julia
# Extract dimers from the full simulation
dimer_systems = get_dimers(systems)

# Generate images showing only dimers
dimer_images = gen_dimer_images(systems, psf;
    photons=1500.0,
    frame_integration=5
)
```

### Analyzing Dimer Formation

To quantify the formation of dimers over time:

```julia
# Calculate fraction of molecules in dimer state
dimer_fractions = analyze_dimer_fraction(systems)

# Plot dimer formation over time
using CairoMakie
times = range(0, params.t_max, length=length(systems))

fig = Figure()
ax = Axis(fig[1, 1],
    xlabel="Time (s)",
    ylabel="Fraction of molecules in dimers",
    title="Dimer formation dynamics"
)
lines!(ax, times, dimer_fractions)
fig
```

## Advanced Usage

### Extending the Model

The interaction-diffusion framework can be extended for more complex simulations:

1. **Multi-state systems**:
   - Modify the state system to include more than just monomers and dimers
   - Add additional reaction types

2. **Spatial heterogeneity**:
   - Introduce position-dependent reaction rates
   - Add fixed obstacles or binding sites

3. **Custom dynamics**:
   - Implement directed motion (e.g., active transport)
   - Add external forces

### Integration with Other Modules

The diffusion simulation integrates well with other SMLMSim components:

```julia
# Generate super-resolution data from diffusion simulation
camera = IdealCamera(128, 128, 0.1)

# Extract emitters from a simulation time point
emitters = [mol.emitter for mol in systems[50].molecules]

# Create SMLD structure
smld = BasicSMLD(emitters, camera, 1, 1, Dict("source" => "diffusion"))

# Apply photophysics and localization uncertainty
using SMLMSim, MicroscopePSFs

# Define fluorophore using positional constructor
fluor = GenericFluor(1e4, [-10.0 10.0; 1.0 -1.0]) # γ=1e4, k_off=10, k_on=1

# ... rest of the example ...
```

## Examples

### Protein Association Kinetics

This example simulates protein association in a membrane:

```julia
# Membrane protein association simulation
params = SmoluchowskiParams(
    density = 0.3,        # sparse proteins
    diff_monomer = 0.05,  # slow membrane diffusion
    k_off = 0.01,         # stable complexes
    t_max = 60.0          # observe over 1 minute
)

systems = simulate(params)
dimer_fractions = analyze_dimer_fraction(systems)

# Calculate time to 50% dimerization
times = range(0, params.t_max, length=length(systems))
t_half = times[findfirst(df -> df >= 0.5, dimer_fractions)]
```

### Single Particle Tracking

This example generates data for single particle tracking analysis:

```julia
# Fast acquisition for tracking
params = SmoluchowskiParams(
    density = 0.05,       # very sparse for tracking
    diff_monomer = 0.2,   # moderate diffusion
    dt = 0.001,           # high temporal resolution
    t_max = 5.0           # 5 second acquisition
)

systems = simulate(params)

# Generate tracking movie with high framerate
psf = MicroscopePSFs.Gaussian2D(0.15)
camera = IdealCamera(128, 128, 0.1)

tracking_movie = gen_image_sequence(psf, systems;
    photons=2000.0,
    frame_integration=10,
    poisson_noise=true
)
```

## API Reference

For a complete list of functions and types, please see the API Reference section.

```@docs
SmoluchowskiParams
DiffusingMolecule
DiffusingMoleculeSystem
```