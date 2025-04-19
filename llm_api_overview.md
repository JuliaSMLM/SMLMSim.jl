# SMLMSim.jl LLM API Overview

This document provides an overview of the SMLMSim.jl package for LLM-assisted usage, including key types, functions, and common workflows for simulating Single Molecule Localization Microscopy (SMLM) data.

## Package Purpose

SMLMSim provides tools for generating realistic SMLM data with:
- Customizable spatial patterns in 2D and 3D
- Stochastic photophysics (blinking kinetics)
- Realistic localization uncertainty
- Molecular diffusion and interactions
- Microscope image generation

## Package Structure and Module Hierarchy

```
SMLMSim.jl
├── Core/
│   ├── CTMC (Continuous Time Markov Chain)
│   ├── Molecules (photophysical models)
│   ├── Patterns (spatial arrangements)
│   └── Photophysics (blinking kinetics)
├── Static/
│   ├── Parameters
│   ├── Coordinate noise
│   └── Simulation
├── InteractionDiffusion/
│   ├── Types (diffusing emitters)
│   ├── Smoluchowski dynamics
│   ├── Analysis
│   └── Helpers
└── CameraImages/
    ├── Image generation
    └── Noise models
```

## Type Definitions

### Abstract Types

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
```

### Parameter Types

```julia
StaticSMLMParams <: SMLMSimParams
    density::Float64            # Patterns per μm²
    σ_psf::Float64              # PSF width in μm
    minphotons::Int             # Minimum photons for detection
    ndatasets::Int              # Number of datasets
    nframes::Int                # Frames per dataset
    framerate::Float64          # Frames per second
    ndims::Int                  # Dimensions (2 or 3)
    zrange::Vector{Float64}     # Axial range for 3D [min_z, max_z]

DiffusionSMLMParams <: SMLMSimParams  # Sometimes referred to as SmoluchowskiParams in docs
    density::Float64            # Molecules per μm²
    box_size::Float64           # Simulation box size (μm)
    diff_monomer::Float64       # Monomer diffusion coefficient (μm²/s)
    diff_dimer::Float64         # Dimer diffusion coefficient (μm²/s)
    diff_dimer_rot::Float64     # Dimer rotational diffusion (rad²/s)
    k_off::Float64              # Dissociation rate (s⁻¹)
    r_react::Float64            # Reaction radius (μm)
    d_dimer::Float64            # Dimer separation (μm)
    dt::Float64                 # Time step (s)
    t_max::Float64              # Total simulation time (s)
    ndims::Int                  # Dimensions (2 or 3)
    boundary::String            # Boundary condition ("periodic"/"reflecting")
    camera_framerate::Float64   # Frames per second (Hz)
    camera_exposure::Float64    # Exposure time (s)
```

### Pattern Definitions

```julia
Nmer2D <: Pattern2D
    n::Int                      # Number of molecules 
    d::Float64                  # Diameter in μm
    x::Vector{Float64}          # X coordinates
    y::Vector{Float64}          # Y coordinates

Line2D <: Pattern2D
    n::Int                      # Number of molecules
    x::Vector{Float64}          # X coordinates
    y::Vector{Float64}          # Y coordinates
    λ::Float64                  # Linear density (molecules/μm)
    endpoints::Vector{Tuple{Float64,Float64}}  # Start and end points

Nmer3D <: Pattern3D
    n::Int                      # Number of molecules
    d::Float64                  # Diameter in μm
    x::Vector{Float64}          # X coordinates
    y::Vector{Float64}          # Y coordinates
    z::Vector{Float64}          # Z coordinates

Line3D <: Pattern3D
    n::Int                      # Number of molecules
    x::Vector{Float64}          # X coordinates
    y::Vector{Float64}          # Y coordinates
    z::Vector{Float64}          # Z coordinates
    λ::Float64                  # Linear density (molecules/μm)
    endpoints::Vector{Tuple{Float64,Float64,Float64}}  # Start and end points
```

### Molecules

```julia
GenericFluor <: Molecule
    γ::AbstractFloat            # Photon emission rate in Hz
    q::Array{<:AbstractFloat}   # State transition rate matrix
```

### Emitter Types

```julia
# Re-exported from SMLMData
AbstractEmitter                  # Base type for all emitters
Emitter2D{T} <: AbstractEmitter  # 2D emitter with (x,y)
Emitter3D{T} <: AbstractEmitter  # 3D emitter with (x,y,z)
Emitter2DFit{T} <: AbstractEmitter # 2D emitter with uncertainties
Emitter3DFit{T} <: AbstractEmitter # 3D emitter with uncertainties

# Defined in SMLMSim
AbstractDiffusingEmitter <: AbstractEmitter  # Base diffusing emitter

DiffusingEmitter2D{T} <: AbstractDiffusingEmitter
    x::T                         # X position (μm)
    y::T                         # Y position (μm)
    photons::T                   # Photon count
    timestamp::T                 # Simulation time (s)
    frame::Int                   # Camera frame
    dataset::Int                 # Dataset number
    id::Int                      # Molecule ID
    state::Symbol                # :monomer or :dimer
    partner_id::Union{Int,Nothing} # Linked molecule ID for dimers

DiffusingEmitter3D{T} <: AbstractDiffusingEmitter
    x::T                         # X position (μm)
    y::T                         # Y position (μm)
    z::T                         # Z position (μm)
    photons::T                   # Photon count
    timestamp::T                 # Simulation time (s)
    frame::Int                   # Camera frame
    dataset::Int                 # Dataset number
    id::Int                      # Molecule ID
    state::Symbol                # :monomer or :dimer
    partner_id::Union{Int,Nothing} # Linked molecule ID for dimers
```

## Function Signatures

### Core Simulation Functions

```julia
# Main simulation interface
simulate(params::StaticSMLMParams; 
         pattern::Union{Pattern,Nothing}=nothing,
         molecule::Molecule=GenericFluor(; q=[-50.0 50.0; 1e-2 -1e-2]),
         camera::AbstractCamera=IdealCamera(128, 128, 0.1)
) -> Tuple{BasicSMLD, BasicSMLD, BasicSMLD}

simulate(params::DiffusionSMLMParams; 
         photons::Float64=1000.0
) -> Array  # Returns a collection of system states, not a single SMLD

# Pattern generation
uniform2D(ρ::T, p::Pattern2D, field_x::T, field_y::T) where T <: AbstractFloat
    -> Tuple{Vector{T}, Vector{T}}

uniform3D(ρ::T, p::Pattern3D, field_x::T, field_y::T; 
         zrange::Vector{T}=T[-1.0, 1.0]) where T <: AbstractFloat
    -> Tuple{Vector{T}, Vector{T}, Vector{T}}

rotate!(p::Pattern2D, θ::T) where T <: AbstractFloat -> nothing
rotate!(p::Pattern3D, R::Matrix{T}) where T <: AbstractFloat -> nothing
rotate!(p::Pattern3D, α::T, β::T, γ::T) where T <: AbstractFloat -> nothing

# Photophysics functions
intensity_trace(f::GenericFluor, nframes::Int, framerate::Real; state1=1) 
    -> Vector{Float64}

kinetic_model(smld::BasicSMLD, f::Molecule, nframes::Int, framerate::Real;
             ndatasets::Int=1, minphotons=50.0, state1::Union{Int, Symbol}=:equilibrium) 
    -> BasicSMLD

# Localization uncertainty
apply_noise(smld::BasicSMLD, σ_psf::AbstractFloat) -> BasicSMLD  # 2D
apply_noise(smld::BasicSMLD, σ_psf::Vector{<:AbstractFloat}) -> BasicSMLD  # 3D

# Diffusion analysis
get_dimers(systems::Array) -> Array
get_monomers(systems::Array) -> Array
analyze_dimer_fraction(systems::Array) -> Tuple{Vector{Int}, Vector{Float64}}
analyze_dimer_lifetime(systems::Array) -> Float64
track_state_changes(systems::Array) -> Dict{Int, Vector{Tuple{Int, Symbol}}}

# Microscope image generation for static SMLD
gen_images(smld::SMLD, psf::AbstractPSF; 
          dataset::Int=1,
          frames=nothing,
          frame_integration::Int=1, 
          support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}}=Inf,
          sampling::Int=2,
          threaded::Bool=true,
          bg::Float64=0.0,
          poisson_noise::Bool=false,
          camera_noise::Bool=false) 
    -> Array{T, 3} where T<:Real

gen_image(smld::SMLD, psf::AbstractPSF, frame::Int; kwargs...) 
    -> Matrix{T} where T<:Real

# Microscope image generation for diffusion systems
gen_image_sequence(
    psf::AbstractPSF, 
    systems::Array;
    photons::Float64=1000.0,
    bg::Float64=0.0,
    frame_integration::Int=1,
    poisson_noise::Bool=false
) -> Array{T, 3} where T<:Real
```

## Common Workflows

### 1. Basic Static SMLM Simulation

```julia
using SMLMSim

# Create camera with physical pixel size
camera = IdealCamera(128, 128, 0.1)  # 128×128 pixels, 100nm pixels

# Create simulation parameters
params = StaticSMLMParams(
    ρ=1.0,                # patterns per μm²
    σ_psf=0.13,           # PSF width in μm (130nm)
    minphotons=50,        # minimum photons for detection
    ndatasets=10,         # number of independent datasets
    nframes=1000,         # frames per dataset
    framerate=50.0        # frames per second
)

# Run simulation
smld_true, smld_model, smld_noisy = simulate(params, camera=camera)

# Extract coordinates from emitters
x_noisy = [e.x for e in smld_noisy.emitters]
y_noisy = [e.y for e in smld_noisy.emitters]
photons = [e.photons for e in smld_noisy.emitters]

# Group emitters by track_id (same original position)
emitters_by_position = Dict()
for e in smld_noisy.emitters
    if !haskey(emitters_by_position, e.track_id)
        emitters_by_position[e.track_id] = []
    end
    push!(emitters_by_position[e.track_id], e)
end
```

### 2. Custom Pattern Simulation

```julia
# Create a custom pattern
pattern = Nmer2D(n=6, d=0.2)  # hexamer with 200nm diameter

# Rotate the pattern if desired
rotate!(pattern, π/4)  # Rotate 45 degrees

# Run simulation with custom pattern
smld_true, smld_model, smld_noisy = simulate(
    StaticSMLMParams(ρ=0.5),  # 0.5 patterns per μm²
    pattern=pattern,
    camera=IdealCamera(128, 128, 0.1)
)

# Create a custom fluorophore model
fluor = GenericFluor(
    γ=20000.0,                # 20,000 photons/s
    q=[0 10 0.1; 1e-2 0 0; 0 0 0]  # ON, OFF, BLEACHED states
)

# Simulate with custom fluorophore
smld_true, smld_model, smld_noisy = simulate(
    StaticSMLMParams(ρ=0.5),
    pattern=pattern,
    molecule=fluor,
    camera=IdealCamera(128, 128, 0.1)
)
```

### 3. 3D Simulation

```julia
# Create 3D simulation parameters
params = StaticSMLMParams(
    ρ=0.5,                # patterns per μm²
    ndims=3,              # 3D simulation
    zrange=[-2.0, 2.0]    # 4μm axial range
)

# Run simulation with 3D pattern
pattern = Nmer3D(n=8, d=0.3)
smld_true, smld_model, smld_noisy = simulate(
    params,
    pattern=pattern,
    camera=IdealCamera(128, 128, 0.1)
)

# Extract 3D coordinates
x = [e.x for e in smld_noisy.emitters]
y = [e.y for e in smld_noisy.emitters]
z = [e.z for e in smld_noisy.emitters]
```

### 4. Diffusion Simulation

```julia
# Set up parameters for Smoluchowski diffusion simulation
params = DiffusionSMLMParams(
    density = 0.5,            # molecules per μm²
    box_size = 10.0,          # μm
    diff_monomer = 0.1,       # μm²/s
    diff_dimer = 0.05,        # μm²/s
    k_off = 0.2,              # s⁻¹
    r_react = 0.01,           # μm
    d_dimer = 0.05,           # μm
    dt = 0.01,                # s
    t_max = 10.0,             # s
    camera_framerate = 20.0,  # fps
    camera_exposure = 0.04    # s
)

# Run simulation - returns a collection of system states
systems = simulate(params)

# Extract dimers
dimer_systems = get_dimers(systems)

# Extract monomers
monomer_systems = get_monomers(systems)

# Analyze dimer formation
frames, fractions = analyze_dimer_fraction(systems)

# Calculate average dimer lifetime
lifetime = analyze_dimer_lifetime(systems)

# Track state changes of individual molecules
state_history = track_state_changes(systems)
```

### 5. Generating Microscope Images

```julia
using MicroscopePSFs

# For static SMLM data:
# Create a PSF model
psf = MicroscopePSFs.GaussianPSF(0.15)  # 150nm PSF width

# Generate a single image for a specific frame
image = gen_image(smld_noisy, psf, 1; 
    bg=5.0,
    poisson_noise=true
)

# Generate a full image sequence
images = gen_images(smld_noisy, psf; 
    frame_integration=10,
    bg=5.0,
    poisson_noise=true
)

# For diffusion simulation:
# Generate images from diffusion systems
diffusion_images = gen_image_sequence(
    psf, 
    systems,
    photons=1000.0,  # Only for gen_image_sequence
    bg=5.0, 
    frame_integration=10,
    poisson_noise=true
)

# Generate images from only dimers
dimer_images = gen_image_sequence(
    psf, 
    dimer_systems, 
    photons=1000.0,
    bg=5.0,
    frame_integration=10
)
```

## Tips for Working with SMLMSim

1. **Units**: All spatial coordinates are in microns (μm), time in seconds (s), and rates in per second (s⁻¹).

2. **Emitter Tracking**: Use `track_id` to track emitters across frames that originated from the same true position.

3. **Camera Settings**: Create a camera with appropriate pixel dimensions and pixel size to match experimental setups:
   ```julia
   # Square camera with same pixel size in both dimensions
   camera = IdealCamera(512, 512, 0.1)  # 512×512 camera, 100nm pixels
   
   # Rectangular camera with different pixel sizes
   camera = IdealCamera(256, 128, (0.1, 0.1))
   ```

4. **Pattern Customization**: Create custom patterns by inheriting from `Pattern2D` or `Pattern3D`:
   ```julia
   # Example of custom 2D pattern
   mutable struct Grid2D <: Pattern2D
       nx::Int      # columns
       ny::Int      # rows
       dx::Float64  # spacing
       dy::Float64  # spacing
       x::Vector{Float64}
       y::Vector{Float64}
   end
   
   function Grid2D(; nx=3, ny=3, dx=0.1, dy=0.1)
       n = nx * ny
       x = zeros(n)
       y = zeros(n)
       
       idx = 1
       for i in 1:nx, j in 1:ny
           x[idx] = (i - (nx+1)/2) * dx
           y[idx] = (j - (ny+1)/2) * dy
           idx += 1
       end
       
       return Grid2D(nx, ny, dx, dy, x, y)
   end
   ```

5. **Diffusion Parameters**: For realistic biological simulations, typical values are:
   - Small proteins in cytoplasm: 0.1-1 μm²/s
   - Membrane proteins: 0.01-0.1 μm²/s
   - Reaction radius: 0.001-0.02 μm (1-20 nm)
   - Dissociation rates (k_off): 0.001-10 s⁻¹ depending on complex stability

6. **Photophysics**: Customize fluorophore blinking behavior using the rate matrix in `GenericFluor`:
   - Two-state model: `q = [0 k_off; k_on 0]`
   - Three-state model with bleaching: `q = [0 k_off k_bleach; k_on 0 0; 0 0 0]`
   - Duty cycle (fraction of time in ON state): `k_on / (k_on + k_off)`

7. **Metadata**: Each SMLD object contains metadata about the simulation parameters:
   ```julia
   # Access pattern type
   pattern_type = smld_noisy.metadata["pattern_type"]
   
   # Access PSF width
   psf_width = smld_noisy.metadata["psf_width"]
   
   # Access all simulation parameters
   params = smld_noisy.metadata["simulation_parameters"]
   ```

8. **Performance Tips**: For large simulations:
   - Use smaller field of view when possible
   - Process data incrementally by simulating one dataset at a time
   - Optimize photophysics parameters to reduce active emitters per frame
   - For 3D simulations, use appropriate z-range for your study
   - Enable multithreading in Julia: `julia --threads=auto`

9. **PSF Types**: Make sure to use the correct PSF type:
   - For static SMLM simulations: `MicroscopePSFs.GaussianPSF(0.15)`
   - The PSF module from MicroscopePSFs.jl must be imported separately