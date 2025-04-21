# SMLMSim.jl API Overview

This document provides a balanced overview of the SMLMSim.jl package for simulating Single Molecule Localization Microscopy (SMLM) data.

## Type Hierarchy

```julia
AbstractSim                      # Base type for all simulation types
└── SMLMSimParams                # Base type for simulation parameters
    ├── StaticSMLMParams         # Parameters for static simulation
    └── DiffusionSMLMParams      # Parameters for diffusion simulation

Pattern                          # Base type for all patterns
├── Pattern2D                    # Base type for 2D patterns
│   ├── Nmer2D                   # N molecules arranged in a circle
│   └── Line2D                   # Linear arrangement of molecules
└── Pattern3D                    # Base type for 3D patterns
    ├── Nmer3D                   # N molecules in a circle at z=0
    └── Line3D                   # 3D linear arrangement

Molecule                         # Base type for molecule models
└── GenericFluor                 # Generic fluorophore with photophysics

AbstractEmitter                  # Base type for all emitters
├── Emitter2D/3D                 # Basic emitters without uncertainties
├── Emitter2DFit/3DFit           # Emitters with uncertainties
└── AbstractDiffusingEmitter     # Base for diffusing emitters
    ├── DiffusingEmitter2D       # 2D diffusing emitter with track_id and partner_id
    └── DiffusingEmitter3D       # 3D diffusing emitter with track_id and partner_id
```

## Core Parameter Types

```julia
StaticSMLMParams(
    density = 1.0,        # Patterns per μm²
    σ_psf = 0.13,         # PSF width in μm
    minphotons = 50,      # Minimum photons for detection
    nframes = 1000,       # Frames per dataset
    framerate = 50.0,     # Frames per second
    ndims = 2,            # Dimensions (2 or 3)
    zrange = [-1.0, 1.0]  # Axial range for 3D [min_z, max_z]
)

DiffusionSMLMParams(
    density = 0.5,        # Molecules per μm²
    box_size = 10.0,      # μm
    diff_monomer = 0.1,   # μm²/s
    diff_dimer = 0.05,    # μm²/s
    k_off = 0.2,          # s⁻¹
    r_react = 0.01,       # μm reaction radius
    d_dimer = 0.05,       # μm dimer separation
    dt = 0.01,            # s time step
    t_max = 10.0          # s total simulation time
)
```

## Key Functions

### Simulation

```julia
# Static SMLM simulation
simulate(params::StaticSMLMParams; 
         starting_conditions=nothing,  # Optional custom emitter configuration
         pattern=Nmer2D(), 
         molecule=GenericFluor(1e4, [-10.0 10.0; 0.5 -0.5]), 
         camera=IdealCamera(128, 128, 0.1))
  -> Tuple{BasicSMLD, BasicSMLD, BasicSMLD}  # (true, model, noisy)

# Diffusion simulation
simulate(params::DiffusionSMLMParams; 
         starting_conditions=nothing,  # Optional custom emitter configuration
         photons=1000.0,
         override_count=nothing)       # Optional exact number of molecules
  -> BasicSMLD  # Single SMLD with all diffusing emitters
```

### Pattern Generation

```julia
# Create patterns
nmer = Nmer2D(n=6, d=0.2)  # 6 molecules in 200nm circle
line = Line2D(λ=5.0, endpoints=[(-1.0, 0.0), (1.0, 0.0)])  # 5 mols/μm

# Generate random pattern distributions
x, y = uniform2D(1.0, pattern, 10.0, 10.0)  # 1 pattern/μm² in 10×10μm field
x, y, z = uniform3D(1.0, pattern3d, 10.0, 10.0, zrange=[-1.0, 1.0])

# Rotate patterns
rotate!(pattern, π/4)  # Rotate 2D pattern by 45 degrees
rotate!(pattern3d, α, β, γ)  # Rotate 3D pattern with Euler angles
```

### Starting Conditions Helpers

```julia
# Convert static emitters to diffusing emitters
diffusing_emitters = convert_to_diffusing_emitters(static_emitters, 1000.0, :monomer)

# Extract final state from a simulation for continued simulation
final_state = extract_final_state(smld)

# Use final state as starting conditions
smld_continued = simulate(params; starting_conditions=final_state)
```

### Photophysics

```julia
# Create fluorophore models
fluor = GenericFluor(1e4, [-10.0 10.0; 0.5 -0.5])  # Positional constructor with γ and q matrix
fluor = GenericFluor(photons=1e4, k_off=10.0, k_on=0.5)  # Keyword constructor for 2-state model

# Generate photophysical data
trace = intensity_trace(fluor, 1000, 50.0)  # Intensity trace for 1000 frames at 50fps
smld_model = kinetic_model(smld_true, fluor, 1000, 50.0)  # Apply blinking model
```

### Localization Uncertainty

```julia
# Add position noise
smld_noisy = apply_noise(smld_model, 0.13)  # 2D: 130nm PSF width
smld_noisy = apply_noise(smld_model, [0.13, 0.13, 0.39])  # 3D: [σx, σy, σz] in μm
```

### Track Utilities

```julia
track_smld = get_track(smld, 5)  # Get emitters from track ID 5
n_tracks = get_num_tracks(smld)  # Count unique tracks
all_tracks = get_tracks(smld)  # Split SMLD by track_id
```

### Diffusion Analysis

```julia
dimer_smld = get_dimers(smld)  # Extract only dimer emitters
monomer_smld = get_monomers(smld)  # Extract only monomer emitters
frames, fractions = analyze_dimer_fraction(smld)  # Dimer fraction per frame
avg_lifetime = analyze_dimer_lifetime(smld)  # Average dimer lifetime
state_history = track_state_changes(smld)  # Track state transitions
```

### Image Generation

```julia
# Generate image stack
images = gen_images(smld, psf;
    frame_integration = 1,    # Time steps per frame
    bg = 0.0,                 # Background photons
    poisson_noise = false     # Add shot noise
)

# Generate single frame
frame_image = gen_image(smld, psf, 10)  # Image of frame 10
```

## Workflow Explanations

### Static SMLM Simulation Workflow

1. **Define simulation parameters**: Create `StaticSMLMParams` with density, PSF width, frames, etc.
2. **Define a molecular pattern**: Create a pattern like `Nmer2D` or `Line2D`
3. **Define a fluorophore model**: Create a `GenericFluor` with emission rate and state transition matrix
4. **Run simulation**: Call `simulate()` to get true, model, and noisy positions
5. **Generate images**: Use `gen_images()` with a PSF model to create synthetic microscope images

```julia
# Complete static simulation example
params = StaticSMLMParams(density=1.0, σ_psf=0.13, nframes=1000)
camera = IdealCamera(128, 128, 0.1)  # 128×128 pixels, 100nm pixels
pattern = Nmer2D(n=6, d=0.2)  # Hexamer with 200nm diameter
molecule = GenericFluor(1e4, [-50.0 50.0; 1e-2 -1e-2])  # Blinking model

# Run simulation
smld_true, smld_model, smld_noisy = simulate(params, pattern=pattern, molecule=molecule, camera=camera)

# Generate images (requires MicroscopePSFs)
using MicroscopePSFs
psf = GaussianPSF(0.15)  # 150nm PSF width
images = gen_images(smld_model, psf, bg=5.0, poisson_noise=true)
```

### Starting with Custom Initial Positions

```julia
# Create custom initial positions
custom_emitters = [
    Emitter2DFit{Float64}(
        x, y,                # Custom positions
        1.0, 0.0,            # photons, background
        0.0, 0.0, 0.0, 0.0;  # uncertainties
        track_id=i           # track ID
    ) for (i, (x, y)) in enumerate(zip([0.0, 0.5], [0.0, 0.0]))
]

# Run simulation with custom initial positions
params = StaticSMLMParams()
smld_true, smld_model, smld_noisy = simulate(
    params;
    starting_conditions=custom_emitters
)
```

### Diffusion-Interaction Simulation Workflow

1. **Define diffusion parameters**: Create `DiffusionSMLMParams` with diffusion coefficients, box size, etc.
2. **Run simulation**: Call `simulate()` to generate diffusing emitters with dimerization
3. **Analyze results**: Extract dimers/monomers and analyze dynamics
4. **Generate images**: Create time-lapse microscopy images with integrated frames

```julia
# Complete diffusion simulation example
params = DiffusionSMLMParams(
    density = 0.5,           # molecules per μm²
    box_size = 10.0,         # μm box size
    diff_monomer = 0.1,      # μm²/s monomer diffusion
    diff_dimer = 0.05,       # μm²/s dimer diffusion
    k_off = 0.2,             # s⁻¹ dissociation rate
    r_react = 0.01,          # μm reaction radius
    d_dimer = 0.05,          # μm dimer separation
    dt = 0.01,               # s simulation time step
    t_max = 10.0,            # s total simulation time
    camera_framerate = 10.0  # fps camera frame rate
)

# Run simulation
smld = simulate(params)

# Analyze dimer formation
frames, dimer_fractions = analyze_dimer_fraction(smld)
dimer_lifetime = analyze_dimer_lifetime(smld)

# Generate time-lapse images
psf = GaussianPSF(0.15)  # 150nm PSF width
images = gen_images(smld, psf, frame_integration=10, bg=5.0, poisson_noise=true)
```

### Multi-Stage Diffusion Simulation with Starting Conditions

```julia
# Initial simulation
params_initial = DiffusionSMLMParams(
    diff_monomer = 0.1,      # μm²/s
    t_max = 5.0              # 5 seconds simulation
)
smld_initial = simulate(params_initial)

# Extract final state
final_state = extract_final_state(smld_initial)

# Continue simulation with different parameters
params_continued = DiffusionSMLMParams(
    diff_monomer = 0.2,      # μm²/s (faster diffusion)
    t_max = 10.0,            # 10 seconds more simulation
    k_off = 0.1              # s⁻¹ (different dissociation rate)
)
smld_continued = simulate(params_continued; starting_conditions=final_state)
```

Notes:
- All spatial coordinates use microns (μm)
- Time is in seconds (s)
- Rates are in per second (s⁻¹)
- The fluorophore rate matrix `q` follows the structure: `q[i,j]` is the transition rate from state i to j, and diagonal elements `q[i,i]` are negative sums of the row
