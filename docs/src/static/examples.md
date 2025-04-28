```@meta
CurrentModule = SMLMSim
```

# Static SMLM Examples

This page provides complete examples for using the static SMLM simulation capabilities of SMLMSim.

## Basic 2D Simulation

This example demonstrates how to simulate and visualize a basic 2D SMLM dataset.

```@example
using SMLMSim
using CairoMakie

# Create camera with physical pixel size
camera = IdealCamera(128, 64, 0.1)  # 128×64 pixels, 100nm pixels

# Create simulation parameters
params = StaticSMLMParams(
    density = 1.0,        # patterns per μm²
    σ_psf = 0.13,         # PSF width in μm
)

# Run simulation with specified pattern and camera
smld_true, smld_model, smld_noisy = simulate(
    params;
    pattern=Nmer2D(n=6, d=0.2),  # hexamer with 200nm diameter
    camera=camera
)

# Extract coordinates and properties from emitters
x_noisy = [e.x for e in smld_noisy.emitters]
y_noisy = [e.y for e in smld_noisy.emitters]
photons = [e.photons for e in smld_noisy.emitters]
frame_nums = [e.frame for e in smld_noisy.emitters]

# Create figure for visualization
fig = Figure(size=(900, 600))

# Super-resolution scatter plot
ax1 = Axis(fig[1, 1], 
    title="Simulated SMLM Localizations",
    xlabel="x (μm)",
    ylabel="y (μm)",
    aspect=DataAspect(),
    yreversed=true  # Makes (0,0) at top-left
)

# Scatter plot with photon counts as color
scatter!(ax1, x_noisy, y_noisy, 
    color=photons,
    colormap=:viridis,
    markersize=4,
    alpha=0.6
)

# Return the figure
fig

```

## 3D Simulation

This example shows how to create and visualize 3D SMLM data.

```@example
using SMLMSim
using CairoMakie

# Create camera with physical pixel size
camera = IdealCamera(128, 128, 0.1)  # 128×128 pixels, 100nm pixels

# Create 3D simulation parameters
params = StaticSMLMParams(
    density = 0.5,        # emitters per μm²
    σ_psf = 0.13,         # PSF width in μm
    ndims = 3,            # 3D simulation
    zrange = [-1.0, 1.0]  # 2μm axial range
)

# Run 3D simulation
smld_true, smld_model, smld_noisy = simulate(
    params;
    pattern=Nmer3D(n=8, d=0.2),  # 3D pattern with 200nm diameter
    camera=camera
)

# Extract 3D coordinates
x = [e.x for e in smld_noisy.emitters]
y = [e.y for e in smld_noisy.emitters]
z = [e.z for e in smld_noisy.emitters]
photons = [e.photons for e in smld_noisy.emitters]

# Create 3D visualization
fig = Figure(size=(900, 700))

ax3d = Axis3(fig[1, 1], 
    xlabel="x (μm)",
    ylabel="y (μm)",
    zlabel="z (μm)",
    title="3D SMLM Simulation",
    aspect = :data
)

# Plot localizations with z-dependent color
scatter!(ax3d, x, y, z,
    color=z,
    colormap=:viridis,
    markersize=5,
    alpha=0.6
)

fig

```

## Generating Microscope Images

This example shows how to generate synthetic microscope images from SMLM simulations:

```@example
using SMLMSim
using MicroscopePSFs
using CairoMakie

# Create camera with physical pixel size
camera = IdealCamera(128, 128, 0.1)  # 128×128 pixels, 100nm pixels

# Create simulation parameters
params = StaticSMLMParams(
    density = 2.0,        # patterns per μm²
    σ_psf = 0.13,         # 130nm PSF width
    nframes = 100,        # 100 frames
    framerate = 20.0      # 20 fps
)

# Run simulation
smld_true, smld_model, smld_noisy = simulate(
    params;
    pattern=Nmer2D(n=6, d=0.2),  # hexamer with 200nm diameter
    camera=camera
)

# Create a PSF model (Gaussian with 150nm width)
psf = MicroscopePSFs.GaussianPSF(0.15)  # 150nm PSF width

# Note: We use smld_model (not smld_noisy) to avoid double-counting uncertainty
# smld_noisy already contains localization errors, and rendering camera images
# naturally introduces noise, so using smld_noisy would apply noise twice
# Generate image stack from emitter data with Poisson noise
images = gen_images(smld_model, psf;
    bg=5.0,            # background photons per pixel
    poisson_noise=true  # add realistic shot noise
)

# Display a single frame
fig = Figure(size=(800, 400))
ax1 = Axis(fig[1, 1], 
    title="Simulated SMLM Image - Frame 20",
    aspect=DataAspect(),
    yreversed=true
)

# Show the 20th frame with inferno colormap
heatmap!(ax1, transpose(images[:, :, 20]), colormap=:inferno)

# Also display emitters for this frame
ax2 = Axis(fig[1, 2], 
    title="Emitter Positions - Frame 20",
    xlabel="x (μm)",
    ylabel="y (μm)",
    aspect=DataAspect()
)

# Extract emitters in frame 20
frame20_emitters = filter(e -> e.frame == 20, smld_model.emitters)
x = [e.x for e in frame20_emitters]
y = [e.y for e in frame20_emitters]
photons = [e.photons for e in frame20_emitters]

# Plot emitter positions
scatter!(ax2, x, y, 
    color=photons,
    colormap=:viridis,
    markersize=10,
    alpha=0.8
)

# Return the figure
fig

```

The resulting `images` is a 3D array with dimensions `[height, width, frames]` that can be used for visualization, algorithm testing, or benchmarking localization software.

## Advanced: Custom Pattern

This example demonstrates creating a custom pattern type for simulation.

```@example
using SMLMSim
using CairoMakie
using Distributions

# Define a custom pattern type: Grid with random jitter
mutable struct JitteredGrid2D <: Pattern2D
    n::Int       # Total number of molecules
    nx::Int      # Columns
    ny::Int      # Rows
    dx::Float64  # Column spacing
    dy::Float64  # Row spacing
    jitter::Float64  # Random position jitter
    x::Vector{Float64}  # x positions
    y::Vector{Float64}  # y positions
end

function JitteredGrid2D(; nx=5, ny=5, dx=0.05, dy=0.05, jitter=0.01)
    n = nx * ny
    x = zeros(n)
    y = zeros(n)
    
    idx = 1
    for i in 1:nx, j in 1:ny
        # Calculate regular grid position
        grid_x = (i - (nx+1)/2) * dx
        grid_y = (j - (ny+1)/2) * dy
        
        # Add random jitter
        jitter_x = rand(Normal(0, jitter))
        jitter_y = rand(Normal(0, jitter))
        
        # Store jittered position
        x[idx] = grid_x + jitter_x
        y[idx] = grid_y + jitter_y
        idx += 1
    end
    
    return JitteredGrid2D(n, nx, ny, dx, dy, jitter, x, y)
end

# Create camera
camera = IdealCamera(64, 64, 0.1)

# Create custom pattern
grid = JitteredGrid2D()

# Create simulation parameters
params = StaticSMLMParams(
    density = 0.2,        # Patterns per μm²
    nframes = 1000        # Number of frames
)

# Run simulation with custom pattern
smld_true, smld_model, smld_noisy = simulate(
    params;
    pattern=grid,
    camera=camera
)

# Visualization
fig = Figure(size=(800, 600))

# Plot ground truth vs. noisy localizations
ax1 = Axis(fig[1, 1], 
    title="Ground Truth",
    xlabel="x (μm)",
    ylabel="y (μm)",
    aspect=DataAspect()
)

ax2 = Axis(fig[1, 2], 
    title="Noisy Localizations",
    xlabel="x (μm)",
    ylabel="y (μm)",
    aspect=DataAspect()
)

# Extract coordinates
x_true = [e.x for e in smld_true.emitters]
y_true = [e.y for e in smld_true.emitters]

x_noisy = [e.x for e in smld_noisy.emitters]
y_noisy = [e.y for e in smld_noisy.emitters]
photons = [e.photons for e in smld_noisy.emitters]

# Plot
scatter!(ax1, x_true, y_true, color=:black, markersize=6)
scatter!(ax2, x_noisy, y_noisy, color=photons, colormap=:viridis, markersize=3, alpha=0.5)

# Return the figure
fig

```
