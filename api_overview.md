# SMLMSim.jl API Overview

SMLMSim is a Julia package for simulating Single Molecule Localization Microscopy (SMLM) data with realistic physical properties. It provides tools for generating static SMLM simulations, diffusion-interaction simulations, and microscope image generation.

## Key Concepts

### Physical Units
All simulations use consistent physical units:
- Spatial dimensions: microns (μm)
- Time: seconds (s)
- Diffusion coefficients: μm²/s
- Rate constants: s⁻¹

### Core Components
- **Patterns**: Spatial arrangements of molecules (e.g., oligomers, lines)
- **Molecules**: Photophysical models of fluorophores with state transitions
- **Simulation**: Parameters and functions for different simulation types
- **Noise**: Realistic localization uncertainty based on photon statistics
- **Camera Images**: Generation of microscope images from simulation data

## Type Hierarchy

- `AbstractSim`: Base type for all simulation types
  - `SMLMSimParams`: Base type for simulation parameters
    - `StaticSMLMParams`: Parameters for static SMLM simulation
    - `DiffusionSMLMParams`: Parameters for diffusion simulation

- `Pattern`: Base type for all molecular patterns
  - `Pattern2D`: Base type for 2D patterns
    - `Nmer2D`: N molecules arranged in a circle
    - `Line2D`: Molecules arranged along a line
  - `Pattern3D`: Base type for 3D patterns
    - `Nmer3D`: N molecules arranged in a circle in 3D
    - `Line3D`: Molecules arranged along a 3D line

- `AbstractLabeling`: Base type for labeling strategies
  - `FixedLabeling`: Deterministic number of fluorophores per site
  - `PoissonLabeling`: Poisson-distributed number of fluorophores per site
  - `BinomialLabeling`: Binomial-distributed number of fluorophores per site

- `Molecule`: Base type for all photophysical models
  - `GenericFluor`: General fluorophore with kinetic state model

- `AbstractDiffusingEmitter`: Base type for diffusing emitters
  - `DiffusingEmitter2D`: 2D emitter with diffusion state
  - `DiffusingEmitter3D`: 3D emitter with diffusion state

## Essential Types

### StaticSMLMParams

Parameters for static SMLM simulation with fixed molecular patterns.

```julia
Base.@kwdef mutable struct StaticSMLMParams <: SMLMSimParams
    density::Float64 = 1.0          # density in particles per square micron
    σ_psf::Float64 = 0.13           # PSF width in microns
    minphotons::Int = 50            # minimum photons for detection
    ndatasets::Int = 10             # number of datasets to simulate
    nframes::Int = 1000             # number of frames per dataset
    framerate::Float64 = 50.0       # frames per second
    ndims::Int = 2                  # dimensionality (2 or 3)
    zrange::Vector{Float64} = [-1.0, 1.0]  # axial range for 3D simulations [min_z, max_z]
end
```

### DiffusionSMLMParams

Parameters for diffusion-based SMLM simulation using Smoluchowski dynamics.

```julia
Base.@kwdef mutable struct DiffusionSMLMParams <: SMLMSimParams
    density::Float64 = 1.0          # number density (molecules/μm²)
    box_size::Float64 = 10.0        # simulation box size (μm)
    diff_monomer::Float64 = 0.1     # monomer diffusion coefficient (μm²/s)
    diff_dimer::Float64 = 0.05      # dimer diffusion coefficient (μm²/s)
    diff_dimer_rot::Float64 = 0.5   # dimer rotational diffusion coefficient (rad²/s)
    k_off::Float64 = 0.2            # dimer dissociation rate (s⁻¹)
    r_react::Float64 = 0.01         # reaction radius (μm)
    d_dimer::Float64 = 0.05         # monomer separation in dimer (μm)
    dt::Float64 = 0.01              # time step (s)
    t_max::Float64 = 10.0           # total simulation time (s)
    ndims::Int = 2                  # number of dimensions (2 or 3)
    boundary::String = "periodic"   # boundary condition type ("periodic" or "reflecting")
    camera_framerate::Float64 = 10.0 # camera frames per second (Hz)
    camera_exposure::Float64 = 0.1   # camera exposure time per frame (s)
end
```

### Molecular Patterns

#### Nmer2D

N molecules symmetrically organized around a circle with diameter d.

```julia
mutable struct Nmer2D <: Pattern2D
    n::Int               # Number of molecules in the pattern
    d::Float64           # Diameter of the circle in microns
    x::Vector{Float64}   # X positions of molecules in microns
    y::Vector{Float64}   # Y positions of molecules in microns
end
```

#### Line2D

Points with uniform random distribution between two endpoints.

```julia
mutable struct Line2D <: Pattern2D
    n::Int                # Number of molecules in the pattern
    x::Vector{Float64}    # X positions of molecules in microns
    y::Vector{Float64}    # Y positions of molecules in microns
    λ::Float64            # Linear molecule density (molecules per micron)
    endpoints::Vector{Tuple{Float64,Float64}}  # Vector of endpoint coordinates
end
```

### Fluorophore Models

#### GenericFluor

Defines a fluorophore with photophysical properties.

```julia
struct GenericFluor <: Molecule
    γ::AbstractFloat     # Photon emission rate in Hz
    q::Array{<:AbstractFloat}  # Rate matrix for state transitions
end
```

### Labeling Strategies

Labeling is separate from photophysics (Molecule). These types control how many fluorophores
attach to each binding site, while Molecule handles blinking kinetics of each fluorophore.

#### AbstractLabeling

Abstract type for labeling strategies. All implementations support an `efficiency` parameter
controlling the probability that a binding site gets labeled at all.

#### FixedLabeling

Deterministic labeling with exactly `n` fluorophores per site.

```julia
struct FixedLabeling <: AbstractLabeling
    n::Int               # Number of fluorophores per site
    efficiency::Float64  # Probability that a site gets labeled (0 to 1)
end
```

#### PoissonLabeling

Poisson-distributed number of fluorophores per site.

```julia
struct PoissonLabeling <: AbstractLabeling
    mean::Float64        # Mean number of fluorophores per site (λ for Poisson)
    efficiency::Float64  # Probability that a site gets labeled (0 to 1)
end
```

Note: With Poisson statistics, some sites may receive 0 fluorophores even when
efficiency=1.0 (especially for small mean values).

#### BinomialLabeling

Binomial-distributed number of fluorophores per site. Models scenarios where each
binding site has `n` potential attachment points, each occupied with probability `p`.

```julia
struct BinomialLabeling <: AbstractLabeling
    n::Int               # Number of potential attachment points per site
    p::Float64           # Probability each attachment point is occupied (0 to 1)
    efficiency::Float64  # Probability that a site gets labeled (0 to 1)
end
```

## Constructor Examples

### Creating Simulation Parameters

```julia
# Static SMLM with default parameters
params_static = StaticSMLMParams()

# Custom parameters for static simulation
params_static = StaticSMLMParams(
    density = 2.0,        # 2 patterns per μm²
    σ_psf = 0.15,         # 150nm PSF width
    nframes = 2000,       # 2000 frames
    framerate = 20.0      # 20 fps
)

# Diffusion simulation with default parameters
params_diff = DiffusionSMLMParams()

# Custom parameters for diffusion simulation
params_diff = DiffusionSMLMParams(
    density = 0.5,        # molecules per μm²
    box_size = 15.0,      # 15μm box size
    diff_monomer = 0.2,   # 0.2 μm²/s diffusion coefficient
    k_off = 0.1,          # 0.1 s⁻¹ dissociation rate
    t_max = 20.0          # 20s simulation time
)
```

### Creating Patterns

```julia
# Create an 8-molecule pattern with 100nm diameter (default)
nmer = Nmer2D()

# Create a custom pattern with 6 molecules and 200nm diameter
hexamer = Nmer2D(n=6, d=0.2)  # d is in microns

# Create a 3D pattern
octamer3d = Nmer3D(n=8, d=0.15)

# Create a line pattern
line = Line2D(λ=5.0, endpoints=[(-2.0, 0.0), (2.0, 0.0)])  # 5 molecules/μm

# Create a 3D line pattern
line3d = Line3D(λ=8.0, endpoints=[(-1.0, 0.0, -0.5), (1.0, 0.0, 0.5)])
```

### Creating Fluorophore Models

```julia
# Create a fluorophore with default parameters
fluor = GenericFluor()

# Create a fluorophore with custom parameters
# γ=100,000 photons/s, k_off=10 Hz, k_on=0.1 Hz
fluor = GenericFluor(1e5, [-10.0 10.0; 0.1 -0.1])

# Create a fluorophore using the 2-state keyword constructor
fluor = GenericFluor(photons=1e5, k_off=50.0, k_on=1e-2)
```

### Creating Labeling Strategies

```julia
# Perfect labeling - exactly 1 fluorophore per site (default)
labeling = FixedLabeling()

# 2 fluorophores per site, 90% of sites labeled
labeling = FixedLabeling(2; efficiency=0.9)

# Poisson labeling - average 1.5 fluorophores per site
labeling = PoissonLabeling(1.5)

# Poisson with 80% labeling efficiency
labeling = PoissonLabeling(2.0; efficiency=0.8)

# Binomial labeling - antibody with 4 dye attachment points, 80% occupancy
labeling = BinomialLabeling(4, 0.8)

# Same with 90% of sites receiving an antibody
labeling = BinomialLabeling(4, 0.8; efficiency=0.9)
```

### Creating a Camera

```julia
# IdealCamera: Poisson noise only
camera = IdealCamera(128, 128, 0.1)  # 128×128 pixels, 100nm pixels

# Specify field of view with an array of pixel edges
pixel_edges_x = 0.0:0.1:12.8  # 0 to 12.8μm in 0.1μm steps
pixel_edges_y = 0.0:0.1:12.8
camera = IdealCamera(pixel_edges_x, pixel_edges_y)

# SCMOSCamera: Realistic noise model with per-pixel calibration (SMLMData 0.4+)
# Parameters: width, height, pixel_size, readnoise
camera_scmos = SCMOSCamera(128, 128, 0.1, 1.6)  # 1.6 e⁻ RMS read noise

# Advanced: Specify all calibration parameters
# offset (ADU), gain (e⁻/ADU), readnoise (e⁻ RMS), quantum efficiency (0-1)
camera_scmos = SCMOSCamera(
    128, 128, 0.1;
    offset=100.0,      # 100 ADU dark level
    gain=0.5,          # 0.5 e⁻/ADU
    readnoise=1.6,     # 1.6 e⁻ RMS
    qe=0.95            # 95% quantum efficiency
)
```

## Core Functions

### Simulation

#### simulate

The main simulation function with multiple methods for different simulation types.

```julia
# Static SMLM simulation
# First create simulation parameters
params = StaticSMLMParams()

# Then run simulation
smld_true, smld_model, smld_noisy = simulate(
    params;
    pattern=Nmer2D(),
    labeling=FixedLabeling(),      # Default: 1 fluorophore per site
    molecule=GenericFluor(),
    camera=IdealCamera(128, 128, 0.1)
)

# With Poisson labeling (variable fluorophores per site)
smld_true, smld_model, smld_noisy = simulate(
    params;
    pattern=Nmer2D(n=6, d=0.2),
    labeling=PoissonLabeling(1.5),  # Average 1.5 fluorophores per site
    molecule=GenericFluor()
)

# Diffusion simulation
# First create simulation parameters
params_diff = DiffusionSMLMParams()

# Then run simulation
smld = simulate(
    params_diff;
    photons=1000.0
)
```

#### kinetic_model

Generate kinetic blinking model from existing localization data.

```julia
# Example of how to call kinetic_model
# First create or obtain the required inputs
smld_true = ... # Some BasicSMLD with true positions
fluor = GenericFluor(1e5, [-10.0 10.0; 0.5 -0.5])  # Fluorophore model
nframes = 1000    # Number of frames
framerate = 50.0  # Frames per second

# Call the function
smld_model = kinetic_model(
    smld_true,     # BasicSMLD with emitter positions
    fluor,         # Molecule with kinetic rates
    nframes,       # Number of frames
    framerate;     # Frames per second
    ndatasets=1,   # Number of independent datasets
    minphotons=50.0, # Minimum photons for detection
    state1=:equilibrium  # Initial state sampling
)
```

#### apply_noise

Add localization uncertainty to emitter positions based on photon counts.

```julia
# Example usage of apply_noise

# First, you need smld_model from a previous step
# For example, from the output of kinetic_model()

# For 2D emitters - add noise with 130nm PSF width
smld_noisy = apply_noise(smld_model, 0.13)  # 0.13 μm = 130nm PSF width

# For 3D emitters - specify PSF width in each dimension
smld_noisy_3d = apply_noise(smld_model_3d, [0.13, 0.13, 0.39])  # [x, y, z] widths in μm
```

### Image Generation

#### gen_images

Generate camera images from SMLD data using the specified PSF model.

```julia
images = gen_images(
    smld::SMLD,
    psf::AbstractPSF;
    dataset::Int=1,                # Dataset number to use from SMLD
    frames=nothing,                # Specific frames to generate (default: all frames)
    support=Inf,                   # PSF support region size (see details below)
    sampling=2,                    # Supersampling factor for PSF integration
    threaded=true,                 # Enable multithreading
    bg=0.0,                        # Background signal level (photons per pixel)
    poisson_noise=false,           # Apply Poisson noise only (simple shot noise)
    camera_noise=false             # Apply full camera noise model (requires SCMOSCamera)
                                   # - For SCMOSCamera: QE, Poisson, read noise, gain, offset
                                   # - For IdealCamera: ignored (use poisson_noise instead)
)

# The support parameter controls PSF computation region:
# 1. Inf (default): Compute PSF over entire image (most accurate but slowest)
support=Inf

# 2. Real number: Use circular region with given radius around each emitter
# Typically 3-5× the PSF width is sufficient for accuracy with better performance
support=1.0  # 1.0 µm radius around each emitter

# 3. Tuple (xmin, xmax, ymin, ymax): Explicit region in microns
support=(4.0, 6.0, 4.0, 6.0)  # Only compute PSF within this region

# Example: sCMOS camera with realistic noise
camera_scmos = SCMOSCamera(128, 128, 0.1, 1.6)
smld = BasicSMLD(emitters, camera_scmos, n_frames, n_datasets)
images_scmos = gen_images(smld, psf, bg=10.0, camera_noise=true)
# Applies: QE → Poisson → read noise → gain → offset

# Example: IdealCamera with Poisson noise only
camera_ideal = IdealCamera(128, 128, 0.1)
smld = BasicSMLD(emitters, camera_ideal, n_frames, n_datasets)
images_poisson = gen_images(smld, psf, bg=10.0, poisson_noise=true)
```

#### gen_image

Generate a single frame camera image.

```julia
# Example of generating a single frame image

# First, define variables
smld = ... # Your SMLD data
psf = GaussianPSF(0.15)  # PSF model with 150nm width
frame_number = 10  # The frame you want to generate

# Generate image for a specific frame
single_frame = gen_image(
    smld,          # SMLD data 
    psf,           # PSF model
    frame_number;  # Frame to generate
    support=1.0,   # Same keyword arguments as gen_images
    bg=5.0,
    poisson_noise=true
)
```

### Analysis Functions

#### Diffusion Analysis

```julia
# Example usage of diffusion analysis functions

# First, run a diffusion simulation
params = DiffusionSMLMParams()
smld = simulate(params)

# Extract dimers from diffusion simulation
dimer_smld = get_dimers(smld)

# Extract monomers
monomer_smld = get_monomers(smld)

# Analyze dimer formation over time
frames, fractions = analyze_dimer_fraction(smld)

# Analyze average dimer lifetime
lifetime = analyze_dimer_lifetime(smld)

# Track state changes over time
state_history = track_state_changes(smld)
```

#### Track Utilities

```julia
# Example usage of track utilities

# First, run a simulation
params = StaticSMLMParams()
smld_true, smld_model, smld_noisy = simulate(params)

# Specify a track ID to extract
track_id = 1  # ID of the track to extract

# Get a specific track by ID
track_smld = get_track(smld_noisy, track_id)

# Get number of unique tracks
n_tracks = get_num_tracks(smld_noisy)

# Get all tracks as separate SMLDs
track_smlds = get_tracks(smld_noisy)
```

### Pattern Manipulation

```julia
# Example usage of pattern manipulation

# Create patterns
pattern2d = Nmer2D(n=6, d=0.2)
pattern3d = Nmer3D(n=8, d=0.15)

# Rotate a 2D pattern by 45 degrees
rotate!(pattern2d, π/4)

# Rotate a 3D pattern with Euler angles (ZYZ convention)
rotate!(pattern3d, π/4, π/6, π/3)  # α, β, γ angles

# Rotate a 3D pattern with a rotation matrix
θ = π/2 # 90 degrees
R = [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]  # Z-axis rotation
rotate!(pattern3d, R)
```

### Labeling Functions

#### apply_labeling

Apply labeling strategy to binding site coordinates, expanding each site to the
appropriate number of fluorophore positions.

```julia
# Generate binding site coordinates from pattern distribution
x, y = uniform2D(1.0, Nmer2D(), 10.0, 10.0)  # ~100 binding sites

# Apply Poisson labeling (average 1.5 fluorophores per site)
x_labeled, y_labeled = apply_labeling((x, y), PoissonLabeling(1.5))

# Apply fixed labeling with efficiency
x_labeled, y_labeled = apply_labeling((x, y), FixedLabeling(1; efficiency=0.8))

# Works with 3D coordinates too
x, y, z = uniform3D(1.0, Nmer3D(), 10.0, 10.0)
x_labeled, y_labeled, z_labeled = apply_labeling((x, y, z), BinomialLabeling(4, 0.8))
```

#### n_fluorophores

Sample the number of fluorophores to place at a single binding site. Called internally
by `apply_labeling`, but can be used directly for custom workflows.

```julia
labeling = PoissonLabeling(1.5)

# Sample number of fluorophores for one site
n = n_fluorophores(labeling)  # Returns Int (can be 0)
```

### Pattern Distribution

Generate random pattern distributions in a field.

```julia
# Create patterns
pattern2d = Nmer2D(n=6, d=0.2)
pattern3d = Nmer3D(n=8, d=0.15)

# Generate random pattern distribution in a field
field_x = 10.0 # μm
field_y = 10.0 # μm
density = 1.0  # patterns per μm²

# Get coordinates for 2D distribution
x, y = uniform2D(density, pattern2d, field_x, field_y)

# Get coordinates for 3D distribution
x, y, z = uniform3D(density, pattern3d, field_x, field_y, zrange=[-2.0, 2.0])
```

## Common Workflows

### Static SMLM Simulation Workflow

1. Define simulation parameters
2. Create a pattern (or use default)
3. Choose a labeling strategy (or use default FixedLabeling)
4. Define a fluorophore model (or use default)
5. Run simulation to get true positions, kinetic model, and noisy localizations
6. Generate microscope images or analyze the data

```julia
# 1. Define parameters
params = StaticSMLMParams(
    density = 1.0,
    σ_psf = 0.13,
    nframes = 1000
)

# 2. Create a pattern
pattern = Nmer2D(n=6, d=0.2)  # hexamer with 200nm diameter

# 3. Choose labeling strategy
labeling = PoissonLabeling(1.5)  # Average 1.5 fluorophores per binding site

# 4. Define fluorophore model
fluor = GenericFluor(photons=1e5, k_off=50.0, k_on=1e-2)

# 5. Run simulation
smld_true, smld_model, smld_noisy = simulate(
    params;
    pattern=pattern,
    labeling=labeling,
    molecule=fluor
)

# 6. Create microscope images with efficient PSF support
psf = GaussianPSF(0.15)  # 150nm PSF width
images = gen_images(smld_model, psf;
    support=1.0,         # 1.0 μm radius around each emitter
    poisson_noise=true   # Add realistic photon counting noise
)
```

### Diffusion Simulation Workflow

1. Define diffusion parameters
2. Run simulation to get emitter trajectories
3. Analyze the diffusion and interaction dynamics
4. Generate microscope images

```julia
# 1. Define parameters
params = DiffusionSMLMParams(
    density = 0.5,        # molecules per μm²
    box_size = 10.0,      # μm
    diff_monomer = 0.1,   # μm²/s
    diff_dimer = 0.05,    # μm²/s
    k_off = 0.2,          # s⁻¹
    t_max = 10.0          # s
)

# 2. Run simulation
smld = simulate(params)

# 3. Analyze the results
dimer_smld = get_dimers(smld)
frames, fractions = analyze_dimer_fraction(smld)

# 4. Generate microscope images with realistic sCMOS noise
# Create sCMOS camera matching simulation box
n_pixels = Int(ceil(params.box_size / 0.1))  # 0.1 μm pixels
camera_scmos = SCMOSCamera(n_pixels, n_pixels, 0.1, 1.6)
smld_cam = BasicSMLD(smld.emitters, camera_scmos, smld.n_frames, 1)

psf = GaussianPSF(0.15)  # 150nm PSF width
images = gen_images(smld_cam, psf;
    support=1.0,         # 1.0 μm PSF support radius (faster)
    bg=5.0,              # Background photons per pixel
    camera_noise=true    # Full sCMOS noise model (QE, Poisson, read noise, gain, offset)
)
```

## Complete Examples

### Static SMLM with Image Generation

```julia
using SMLMSim
using MicroscopePSFs

# Define a camera and simulation parameters
camera = IdealCamera(128, 128, 0.1)  # 128×128 pixels, 100nm pixels
params = StaticSMLMParams(density=1.0, σ_psf=0.13, nframes=1000) 

# Run simulation for an 8-molecule ring pattern
smld_true, smld_model, smld_noisy = simulate(
    params;
    pattern=Nmer2D(n=8, d=0.1),  # 100nm diameter ring
    molecule=GenericFluor(1e5, [-10.0 10.0; 0.5 -0.5]),  # γ=100,000, k_off=10, k_on=0.5
    camera=camera
)

# Create a PSF model
psf = GaussianPSF(0.15)  # 150nm PSF width

# Generate microscope images with finite PSF support
images = gen_images(smld_model, psf;
    support=1.0,         # 1.0 μm PSF support radius (faster than Inf)
    bg=5.0,              # 5 background photons per pixel
    poisson_noise=true   # Add realistic photon counting noise
)

println("Generated $(length(smld_noisy.emitters)) localizations and $(size(images,3)) images.")
```

### Diffusion with Dimer Analysis

```julia
using SMLMSim
using MicroscopePSFs

# Set diffusion simulation parameters
params = DiffusionSMLMParams(
    density = 0.5,        # molecules per μm²
    box_size = 10.0,      # μm
    diff_monomer = 0.1,   # μm²/s
    diff_dimer = 0.05,    # μm²/s
    k_off = 0.2,          # s⁻¹
    t_max = 10.0,         # s
    boundary = "reflecting"  # Use reflecting boundaries
)

# Run diffusion simulation
smld = simulate(params)

# Analyze dimer formation
frames, dimer_fraction = analyze_dimer_fraction(smld)
avg_lifetime = analyze_dimer_lifetime(smld)

# Generate microscope images with finite PSF support
psf = GaussianPSF(0.15)  # 150nm PSF width
images = gen_images(smld, psf; 
    support=1.0,         # 1.0 μm PSF support radius (faster)
    bg=2.0,              # Background photons per pixel
    poisson_noise=true   # Add realistic photon counting noise
)

println("Simulation complete with $(length(smld.emitters)) emitters")
println("Average dimer fraction: $(mean(dimer_fraction))")
println("Average dimer lifetime: $(avg_lifetime) seconds")
```

### Custom Pattern with 3D Simulation

```julia
using SMLMSim
using MicroscopePSFs

# Define a custom 3D pattern: two rings at different z-positions
mutable struct DoubleRing3D <: Pattern3D
    n::Int
    d1::Float64
    d2::Float64
    z1::Float64
    z2::Float64
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end

function DoubleRing3D(; n=8, d1=0.1, d2=0.2, z1=-0.2, z2=0.2)
    total_n = 2*n
    x = zeros(total_n)
    y = zeros(total_n)
    z = zeros(total_n)
    
    # First ring (bottom)
    for i = 1:n
        θ = 2π * (i-1) / n
        x[i] = d1/2 * cos(θ)
        y[i] = d1/2 * sin(θ)
        z[i] = z1
    end
    
    # Second ring (top)
    for i = 1:n
        θ = 2π * (i-1) / n + π/n  # Offset angle for second ring
        x[n+i] = d2/2 * cos(θ)
        y[n+i] = d2/2 * sin(θ)
        z[n+i] = z2
    end
    
    return DoubleRing3D(n, d1, d2, z1, z2, x, y, z)
end

# Create camera and parameters
camera = IdealCamera(128, 128, 0.1)
params = StaticSMLMParams(
    density = 0.5,
    σ_psf = 0.13,
    nframes = 2000,
    ndims = 3,
    zrange = [-1.0, 1.0]
)

# Create custom pattern
double_ring = DoubleRing3D(n=6, d1=0.15, d2=0.3, z1=-0.3, z2=0.3)

# Run simulation
smld_true, smld_model, smld_noisy = simulate(
    params;
    pattern=double_ring,
    camera=camera
)

# Generate images with a 3D astigmatic PSF and finite support
# Create a PSF with astigmatism using Zernike coefficients
using MicroscopePSFs
zc = ZernikeCoefficients(15)
zc.phase[6] = 0.5  # Add vertical astigmatism
psf_scalar = ScalarPSF(1.4, 0.532, 1.52; zernike_coeffs=zc)

# Create SplinePSF for speed
xy_sampling, z_sampling = 0.05, 0.1
x_range = y_range = -1.0:xy_sampling:1.0
z_range = -1.0:z_sampling:1.0
psf_spline = SplinePSF(psf_scalar, x_range, y_range, z_range)

# Generate images using the spline PSF with finite support
images = gen_images(smld_model, psf_spline; 
    support=0.5,         # 0.5 μm PSF support radius for performance
    bg=5.0,              # Background photons per pixel
    poisson_noise=true   # Add realistic photon counting noise
)

println("Generated $(length(smld_noisy.emitters)) localizations in 3D")
```