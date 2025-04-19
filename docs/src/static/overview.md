```@meta
CurrentModule = SMLMSim
```

# Static SMLM Simulation Overview

Static SMLM simulation creates fixed molecular patterns with realistic blinking kinetics and localization uncertainty. This approach is ideal for simulating structured biological samples like protein complexes, filaments, and other cellular structures.

## Simulation Workflow

Static simulations follow this general workflow:

1. **Pattern definition**: Define molecular arrangements (e.g., oligomers, lines)
2. **Pattern distribution**: Distribute patterns throughout the field of view
3. **Photophysics**: Apply stochastic blinking behavior to each molecule
4. **Localization uncertainty**: Add realistic position errors based on photon counts

## Simulation Parameters

Static simulations are configured using the `StaticSMLMParams` type:

```julia
# Default parameters
params = StaticSMLMParams()

# Custom parameters
params = StaticSMLMParams(
    ρ=2.0,                # 2 patterns per μm²
    σ_psf=0.15,           # 150nm PSF width
    minphotons=100,       # Minimum photons for detection
    ndatasets=5,          # Number of independent datasets
    nframes=2000,         # Frames per dataset
    framerate=100.0,      # Frames per second
    ndims=3,              # 3D simulation
    zrange=[-2.0, 2.0]    # 4μm axial range
)
```

### Key Parameters

- **`ρ` (density)**: Number of patterns per square micron
- **`σ_psf`**: Point spread function width in microns (impacts localization precision)
- **`minphotons`**: Minimum photons required for detection
- **`nframes`**: Number of frames in the simulation
- **`framerate`**: Frame rate in frames per second
- **`ndims`**: Dimensionality (2 or 3)
- **`zrange`**: Axial range for 3D simulations

## Running a Simulation

The high-level `simulate()` function provides a simple interface for static simulations:

```julia
# Define a camera
camera = IdealCamera(1:128, 1:128, 0.1)  # 128×128 pixels, 100nm pixels

# Run simulation with parameters
smld_true, smld_model, smld_noisy = simulate(
    params,
    pattern=Nmer2D(n=6, d=0.2),  # hexamer with 200nm diameter
    molecule=GenericFluor(γ=1e4, q=[0 10; 0.5 0]),  # fluorophore model
    camera=camera
)
```

Alternatively, you can use keyword arguments directly:

```julia
smld_true, smld_model, smld_noisy = simulate(
    ρ=1.0,                # patterns per μm²
    σ_psf=0.13,           # PSF width in μm
    nframes=1000,         # frames
    framerate=50.0,       # frames per second
    pattern=Nmer2D(n=8, d=0.1),  # pattern type
    camera=camera
)
```

## Understanding Simulation Results

The `simulate()` function returns three `SMLD` objects:

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

Each SMLD object contains emitters with properties like position, frame number, photon count, and uncertainties.

## Working with Results

You can extract information from the simulation results for analysis or visualization:

```julia
# Extract coordinates from noisy emitters
x_noisy = [e.x for e in smld_noisy.emitters]
y_noisy = [e.y for e in smld_noisy.emitters]
photons = [e.photons for e in smld_noisy.emitters]
frames = [e.frame for e in smld_noisy.emitters]

# Group emitters by original position using track_id
emitters_by_position = Dict()
for e in smld_noisy.emitters
    if !haskey(emitters_by_position, e.track_id)
        emitters_by_position[e.track_id] = []
    end
    push!(emitters_by_position[e.track_id], e)
end
```

## 2D vs 3D Simulations

The dimensionality of the simulation is controlled by both the `ndims` parameter and the pattern type:

```julia
# 2D simulation
params = StaticSMLMParams(ndims=2)
pattern2d = Nmer2D(n=6, d=0.2)
smld_true_2d, smld_model_2d, smld_noisy_2d = simulate(params, pattern=pattern2d)

# 3D simulation
params = StaticSMLMParams(ndims=3, zrange=[-2.0, 2.0])
pattern3d = Nmer3D(n=6, d=0.2)
smld_true_3d, smld_model_3d, smld_noisy_3d = simulate(params, pattern=pattern3d)
```

For 3D simulations, the localization uncertainty is automatically scaled in the axial dimension to reflect the typically lower z-resolution in SMLM experiments.

