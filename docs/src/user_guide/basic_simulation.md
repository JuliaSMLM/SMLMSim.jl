```@meta
CurrentModule = SMLMSim
```

# Basic Simulation

This guide covers the fundamental simulation capabilities of SMLMSim, focusing on the high-level `simulate()` function which provides a convenient interface for generating SMLM data.

## Quick Start

The simplest way to generate SMLM data is:

```julia
using SMLMSim

# Create a camera with 100nm pixels
camera = IdealCamera(128, 128, 0.1)

# Run simulation with default parameters
smld_true, smld_model, smld_noisy = simulate(camera=camera)
```

This basic simulation creates:
- 8-molecule circular patterns (Nmer2D with n=8, d=0.1μm)
- 1 pattern per square micron (ρ=1.0)
- PSF width of 130nm (σ_psf=0.13μm)
- Two-state fluorophore with realistic blinking
- 1000 frames at 50 frames per second
- Minimum photon threshold of 50 for detection

## Understanding the Output

The `simulate()` function returns three SMLD objects:

1. **`smld_true`**: Ground truth emitter positions
   - Contains the exact coordinates of all emitters
   - Single frame (no temporal information)
   - No blinking or noise applied

2. **`smld_model`**: Positions with blinking kinetics
   - Subset of true positions appearing in different frames
   - Realistic blinking behavior based on kinetic model
   - No position uncertainty (exact coordinates)

3. **`smld_noisy`**: Positions with blinking and localization uncertainty
   - Same temporal distribution as `smld_model`
   - Position noise based on photon statistics
   - Most comparable to real experimental data

Each SMLD object contains individual emitters with properties like position, frame number, photon count, and uncertainties.

## Customizing Simulation Parameters

For more control, you can customize various parameters:

```julia
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

### Key Parameters

- **`ρ`**: Density of patterns per square micron
- **`σ_psf`**: Point spread function width in microns
- **`minphotons`**: Photon threshold for detection
- **`ndatasets`**: Number of independent datasets to simulate
- **`nframes`**: Number of frames per dataset
- **`framerate`**: Frame rate in frames per second
- **`pattern`**: Pattern type and configuration
- **`molecule`**: Fluorophore model with kinetic parameters
- **`camera`**: Camera configuration
- **`zrange`**: Axial range for 3D simulations

## 2D vs 3D Simulations

SMLMSim can generate both 2D and 3D data:

```julia
# 2D simulation
pattern2d = Nmer2D(n=8, d=0.1)
smld_true_2d, smld_model_2d, smld_noisy_2d = simulate(pattern=pattern2d)

# 3D simulation
pattern3d = Nmer3D(n=8, d=0.1)
smld_true_3d, smld_model_3d, smld_noisy_3d = simulate(
    pattern=pattern3d,
    zrange=[-2.0, 2.0]  # 4μm axial range
)
```

For 3D simulations, the PSF width is automatically adjusted in the axial dimension, typically using a factor of 3× compared to lateral dimensions.

## Working with Simulation Results

You can access individual emitter properties from the SMLD objects:

```julia
# Extract coordinates from noisy emitters
x_noisy = [e.x for e in smld_noisy.emitters]
y_noisy = [e.y for e in smld_noisy.emitters]
photons = [e.photons for e in smld_noisy.emitters]
frames = [e.frame for e in smld_noisy.emitters]
track_ids = [e.track_id for e in smld_noisy.emitters]

# Group emitters by original position using track_id
emitters_by_position = Dict()
for e in smld_noisy.emitters
    if !haskey(emitters_by_position, e.track_id)
        emitters_by_position[e.track_id] = []
    end
    push!(emitters_by_position[e.track_id], e)
end
```

The `track_id` field is particularly useful as it links emitters across frames that originated from the same true position.

## Accessing Simulation Metadata

Each SMLD object contains metadata about the simulation:

```julia
# Get simulation parameters
density = smld_noisy.metadata["density"]
pattern_type = smld_noisy.metadata["pattern_type"]
psf_width = smld_noisy.metadata["psf_width"]
```

This metadata can be useful for keeping track of simulation parameters and for reproducibility.