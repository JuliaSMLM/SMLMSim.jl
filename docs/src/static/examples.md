```@meta
CurrentModule = SMLMSim
```

# Static SMLM Examples

This page provides complete examples for using the static SMLM simulation capabilities of SMLMSim.

## Basic 2D Simulation

This example demonstrates how to simulate and visualize a basic 2D SMLM dataset.

```julia
using SMLMSim
using CairoMakie

# Create camera with physical pixel size
camera = IdealCamera(128, 256, 0.1)  # 128×256 pixels, 100nm pixels

# Create simulation parameters
params = StaticSMLMParams(
    density = 1.0,        # emitters per μm²
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

Colorbar(fig[1, 2], colormap=:viridis, label="Photons")

# Temporal color-coding in second plot
ax2 = Axis(fig[2, 1], 
    title="Temporal Color-Coding",
    xlabel="x (μm)",
    ylabel="y (μm)",
    aspect=DataAspect(),
    yreversed=true
)

scatter!(ax2, x_noisy, y_noisy, 
    color=frame_nums,
    colormap=:plasma,
    markersize=4,
    alpha=0.6
)

Colorbar(fig[2, 2], colormap=:plasma, label="Frame Number")

fig
# output

```

## 3D Simulation

This example shows how to create and visualize 3D SMLM data.

```julia
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
    aspect=(1, 1, 1)
)

# Plot localizations with z-dependent color
scatter!(ax3d, x, y, z,
    color=z,
    colormap=:viridis,
    markersize=5,
    alpha=0.6
)

Colorbar(fig[1, 2], colormap=:viridis, label="z position (μm)")

# Create 2D projections
ax_xy = Axis(fig[2, 1], 
    title="XY Projection",
    xlabel="x (μm)",
    ylabel="y (μm)",
    aspect=DataAspect(),
    yreversed=true
)

ax_xz = Axis(fig[2, 2], 
    title="XZ Projection",
    xlabel="x (μm)",
    ylabel="z (μm)",
    aspect=DataAspect()
)

scatter!(ax_xy, x, y, color=z, colormap=:viridis, markersize=3, alpha=0.6)
scatter!(ax_xz, x, z, color=z, colormap=:viridis, markersize=3, alpha=0.6)

fig
# output

```

## Generating Microscope Images

This example shows how to generate synthetic microscope images from SMLM simulations:

```julia
using SMLMSim
using MicroscopePSFs
using CairoMakie

# Create camera with physical pixel size
camera = IdealCamera(128, 128, 0.1)  # 128×128 pixels, 100nm pixels

# Create simulation parameters
params = StaticSMLMParams(
    density = 0.5,        # patterns per μm²
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

# Generate image stack from emitter data with Poisson noise
images = gen_images(smld_noisy, psf;
    bg=10.0,            # background photons per pixel
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
frame20_emitters = filter(e -> e.frame == 20, smld_noisy.emitters)
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

# Set same limits as the image
limits!(ax2, 0, 12.8, 0, 12.8)  # 128 pixels * 0.1 μm = 12.8 μm

fig
# output

```

The resulting `images` is a 3D array with dimensions `[height, width, frames]` that can be used for visualization, algorithm testing, or benchmarking localization software.

## Customizing Photophysics

This example demonstrates customizing the photophysical properties of fluorophores.

```julia
using SMLMSim
using CairoMakie

# Create camera
camera = IdealCamera(128, 128, 0.1)

# Define different fluorophore models
# 1. Slow blinking (long on/off times)
fluor_slow = GenericFluor(photons=10000.0, k_off=0.5, k_on=0.05)  # Slow rates

# 2. Fast blinking (short on/off times)
fluor_fast = GenericFluor(photons=10000.0, k_off=50.0, k_on=5.0)  # Fast rates

# 3. Three-state model (using direct constructor for advanced models)
# This model has photobleaching as third state
q_bleach = [-10.1 10.0 0.1; 1.0 -1.0 0.0; 0.0 0.0 0.0]  # Third state is absorbing
fluor_bleach = GenericFluor(10000.0, q_bleach)

# Create simulation parameters
params = StaticSMLMParams(
    density = 0.2,
    nframes = 2000
)

# Run simulations with different fluorophore models
_, _, smld_slow = simulate(
    params;
    pattern=Nmer2D(n=3, d=0.2),
    molecule=fluor_slow,
    camera=camera
)

_, _, smld_fast = simulate(
    params;
    pattern=Nmer2D(n=3, d=0.2),
    molecule=fluor_fast,
    camera=camera
)

_, _, smld_bleach = simulate(
    params;
    pattern=Nmer2D(n=3, d=0.2),
    molecule=fluor_bleach,
    camera=camera
)

# Extract data for analysis
function get_time_traces(smld)
    # Group emitters by track_id (original position)
    emitters_by_track = Dict()
    for e in smld.emitters
        if !haskey(emitters_by_track, e.track_id)
            emitters_by_track[e.track_id] = []
        end
        push!(emitters_by_track[e.track_id], e)
    end
    
    # Get a representative trace (first track_id)
    if isempty(emitters_by_track)
        # Return an empty trace if no emitters are found
        return Float64[]
    end
    
    first_track = minimum(keys(emitters_by_track))
    emitters = emitters_by_track[first_track]
    
    # Create frame-by-frame intensity trace
    if isempty(emitters)
        return Float64[]
    end
    
    max_frame = maximum([e.frame for e in emitters])
    trace = zeros(max_frame)
    for e in emitters
        trace[e.frame] = e.photons
    end
    
    return trace
end

# Get traces for each model
trace_slow = get_time_traces(smld_slow)
trace_fast = get_time_traces(smld_fast)
trace_bleach = get_time_traces(smld_bleach)

# Visualization
fig = Figure(size=(900, 700))

# Plot time traces
ax1 = Axis(fig[1, 1], 
    title="Slow Blinking",
    xlabel="Frame",
    ylabel="Photons"
)

ax2 = Axis(fig[2, 1], 
    title="Fast Blinking",
    xlabel="Frame",
    ylabel="Photons"
)

ax3 = Axis(fig[3, 1], 
    title="Photobleaching Model",
    xlabel="Frame",
    ylabel="Photons"
)

# Plot traces
if !isempty(trace_slow)
    stem!(ax1, 1:length(trace_slow), trace_slow)
end
if !isempty(trace_fast)
    stem!(ax2, 1:length(trace_fast), trace_fast)
end
if !isempty(trace_bleach)
    stem!(ax3, 1:length(trace_bleach), trace_bleach)
end

fig
# output

```

## Advanced: Custom Pattern

This example demonstrates creating a custom pattern type for simulation.

```julia
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
camera = IdealCamera(128, 128, 0.1)

# Create custom pattern
grid = JitteredGrid2D(nx=8, ny=8, dx=0.05, dy=0.05, jitter=0.01)

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

fig
# output

```