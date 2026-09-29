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

The main interface for running diffusion simulations is the `simulate` function with `DiffusionSMLMConfig`:

```julia
using SMLMSim

# Set simulation parameters
params = DiffusionSMLMConfig(
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
smld, info = simulate(params)
```

The `smld` output is a `BasicSMLD` structure containing all emitters across all time points, with each emitter having frame information corresponding to the camera settings. The `info` struct contains additional simulation metadata such as `info.elapsed_s`.

## Simulation Parameters

The `DiffusionSMLMConfig` structure allows you to customize various aspects of the simulation:

```julia
# More complex simulation
params = DiffusionSMLMConfig(
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

`dt` is the physics step and also sets the sub-steps per frame (motion blur):
`camera_exposure` and `1/camera_framerate` should be integer multiples of `dt`
(otherwise `simulate` rounds to the nearest step count and warns).
The `γ` argument of `simulate` is the emission rate in photons/s; each of the
`n_sub` records in a frame carries `γ·dt`, so a frame holds `γ·n_sub·dt` photons. That
equals `γ·camera_exposure` when `camera_exposure` is an integer multiple of `dt` and not longer
than the frame period `1/camera_framerate`; otherwise `n_sub = round(camera_exposure/dt)`
(at least 1), capped at the number of steps in the frame period, with a warning. Without `γ`, each record carries 1000 photons (γ = 1000/dt,
0.7's default; 0.8.0 will change the default to a fixed rate). The `photons` keyword is
deprecated (γ = photons/dt) and is removed in 0.8.0.

To continue a run, pass `starting_conditions=smld` or `extract_end_state(smld)`; both resume
at the exact end state of an unchanged run. Brightness: a run whose rate was set with `γ` continues at
that rate (each record carries γ·dt at the new `dt`) when every resumed molecule still carries γ·dt at
the saved `dt`, to a relative 1e-6; otherwise (a default or `photons` source, or an edited or
concatenated SMLD) each molecule keeps its photons per record, as in 0.7.1, with a warning when a γ
rate could not be kept. D: a track keeps its saved D (drawn from a `monomer_mobility` mixture) when its
record is the run's own, with a warning when the new mixture differs from the run's; every other track
draws from the new mixture or uses `diff_monomer` at run time, which is never saved, so a later change
applies. `extract_final_state` is deprecated.

Monomers can be given a mixture of mobility populations with
`monomer_mobility = [(0.85, 0.0), (0.05, 0.08), (0.10, 0.38)]` (entries are
`(fraction, D)`); the drawn coefficient per molecule is stored in
`smld.metadata["monomer_D"]`. `frame_dimer_truth(smld)` returns per-frame, per-molecule
dimer ground truth.


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
psf = MicroscopePSFs.GaussianPSF(0.15)  # 150nm PSF width

# Generate images
image_stack, img_info = gen_images(smld, psf;
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

