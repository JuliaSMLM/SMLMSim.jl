```@meta
CurrentModule = SMLMSim
```

# Examples

This page provides complete workflow examples for using SMLMSim to generate and analyze SMLM data.

## Basic 2D Simulation with Visualization

This example demonstrates how to simulate and visualize a basic 2D SMLM dataset.

```julia
using SMLMSim
using CairoMakie

# Create camera with physical pixel size
camera = IdealCamera(1:128, 1:256, 0.1)  # 128×256 pixels, 100nm pixels

# Simulation parameters in physical units
smld_true, smld_model, smld_noisy = simulate(;
    ρ=1.0,                # emitters per μm²
    σ_psf=0.13,           # PSF width in μm
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

# Display figure
fig
```

## 3D Simulation and Visualization

This example shows how to create and visualize 3D SMLM data.

```julia
using SMLMSim
using CairoMakie

# Create camera with physical pixel size
camera = IdealCamera(1:128, 1:128, 0.1)  # 128×128 pixels, 100nm pixels

# 3D simulation parameters
smld_true, smld_model, smld_noisy = simulate(;
    ρ=0.5,                # emitters per μm²
    σ_psf=0.13,           # PSF width in μm
    pattern=Nmer3D(n=8, d=0.2),  # 3D pattern with 200nm diameter
    camera=camera,
    zrange=[-1.0, 1.0]    # 2μm axial range
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
```

## Analyzing Simulated Data

This example demonstrates how to analyze simulated SMLM data to extract structural information.

```julia
using SMLMSim
using CairoMakie
using Statistics
using StatsBase

# Simulate a hexamer structure
pattern = Nmer2D(n=6, d=0.2)  # 6 molecules in a 200nm circle
camera = IdealCamera(1:128, 1:128, 0.1)
smld_true, smld_model, smld_noisy = simulate(
    pattern=pattern,
    ρ=0.5,       # patterns per μm²
    nframes=500,
    camera=camera
)

# Group localizations by track_id to identify clusters
clusters = Dict()
for e in smld_noisy.emitters
    if !haskey(clusters, e.track_id)
        clusters[e.track_id] = []
    end
    push!(clusters[e.track_id], e)
end

# Function to calculate cluster properties
function analyze_cluster(emitters)
    # Extract coordinates
    x = [e.x for e in emitters]
    y = [e.y for e in emitters]
    
    # Calculate center
    center_x = mean(x)
    center_y = mean(y)
    
    # Calculate distances from center
    distances = [sqrt((x[i] - center_x)^2 + (y[i] - center_y)^2) for i in 1:length(x)]
    
    # Calculate mean diameter and standard deviation
    diameter = 2 * mean(distances)
    diameter_std = 2 * std(distances)
    
    # Count localizations
    n_localizations = length(emitters)
    
    return (
        center_x=center_x,
        center_y=center_y,
        diameter=diameter,
        diameter_std=diameter_std,
        n_localizations=n_localizations
    )
end

# Analyze each cluster
cluster_properties = [analyze_cluster(emitters) for emitters in values(clusters)]

# Filter clusters with enough localizations
min_locs = 10
filtered_clusters = filter(c -> c.n_localizations >= min_locs, cluster_properties)

# Extract diameter statistics
diameters = [c.diameter for c in filtered_clusters]
diameter_stds = [c.diameter_std for c in filtered_clusters]

# Create histogram of measured diameters
fig = Figure(size=(900, 400))

ax1 = Axis(fig[1, 1], 
    title="Measured Cluster Diameters",
    xlabel="Diameter (μm)",
    ylabel="Count"
)

hist!(ax1, diameters, bins=20)

# Mark the true diameter
vlines!(ax1, [pattern.d], color=:red, linestyle=:dash, linewidth=2,
       label="True diameter ($(pattern.d) μm)")

axislegend(ax1)

# Compare measured vs. true diameters
ax2 = Axis(fig[1, 2],
    title="Measurement Precision",
    xlabel="Number of Localizations",
    ylabel="Diameter Error (μm)"
)

n_locs = [c.n_localizations for c in filtered_clusters]
errors = abs.([c.diameter - pattern.d for c in filtered_clusters])

scatter!(ax2, n_locs, errors, markersize=8, alpha=0.6)

fig
```

## Interaction-Diffusion Simulation

This example shows how to simulate molecular interactions using the Smoluchowski dynamics model.

```julia
using SMLMSim
using CairoMakie
using MicroscopePSFs

# Set up parameters for Smoluchowski diffusion simulation
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

# Run diffusion simulation
systems = simulate(params)

# Generate microscope images
# First, set up a camera and PSF
pixelsize = 0.1  # 100nm pixels
pixels = Int64(round(params.box_size/pixelsize))
camera = IdealCamera(; xpixels=pixels, ypixels=pixels, pixelsize=pixelsize)

# Set up PSF (using Airy PSF as an example)
na = 1.3
wavelength = 0.6  # microns
psf = MicroscopePSFs.Gaussian2D(0.15)  # 150nm PSF width

# Generate image sequence with frame integration
images = gen_image_sequence(
    psf, 
    systems,
    photons=1000.0, 
    bg=5.0, 
    frame_integration=10
)

# Extract only dimers
dimer_systems = get_dimers(systems)
dimer_images = gen_image_sequence(
    psf, 
    dimer_systems, 
    photons=1000.0,
    bg=5.0,
    frame_integration=10
)

# Visualize microscope images
function display_microscope_images(images, dimer_images)
    fig = Figure(size=(1000, 400))
    
    # All molecules image
    ax1 = Axis(fig[1, 1], 
        title="All Molecules",
        xlabel="x (μm)",
        ylabel="y (μm)",
        aspect=DataAspect()
    )
    
    # Dimer-only image
    ax2 = Axis(fig[1, 2], 
        title="Dimers Only",
        xlabel="x (μm)",
        ylabel="y (μm)",
        aspect=DataAspect()
    )
    
    # Display example frames
    frame = size(images, 3) ÷ 2  # Middle frame
    heatmap!(ax1, images[:, :, frame], colormap=:inferno)
    heatmap!(ax2, dimer_images[:, :, frame], colormap=:inferno)
    
    return fig
end

display_microscope_images(images, dimer_images)
```

## Advanced: Custom Pattern with Temporal Dynamics

This example demonstrates creating a custom pattern and simulating temporal dynamics.

```julia
using SMLMSim
using CairoMakie

# Create a custom pattern: two converging lines
mutable struct ConvergingLines <: Pattern2D
    n1::Int       # number of molecules in first line
    n2::Int       # number of molecules in second line
    length::Float64  # line length
    angle::Float64   # angle between lines in radians
    x::Vector{Float64}
    y::Vector{Float64}
end

function ConvergingLines(; n1=10, n2=10, length=1.0, angle=π/6)
    n_total = n1 + n2
    x = zeros(n_total)
    y = zeros(n_total)
    
    # Fill first line
    for i in 1:n1
        t = (i - 1) / (n1 - 1)  # parameter from 0 to 1
        x[i] = t * length * cos(-angle/2)
        y[i] = t * length * sin(-angle/2)
    end
    
    # Fill second line
    for i in 1:n2
        t = (i - 1) / (n2 - 1)  # parameter from 0 to 1
        x[n1+i] = t * length * cos(angle/2)
        y[n1+i] = t * length * sin(angle/2)
    end
    
    return ConvergingLines(n1, n2, length, angle, x, y)
end

# Create the custom pattern
pattern = ConvergingLines(n1=15, n2=15, length=1.0, angle=π/3)

# Simulate with custom pattern
camera = IdealCamera(1:128, 1:128, 0.1)
smld_true, smld_model, smld_noisy = simulate(
    pattern=pattern,
    ρ=0.2,       # patterns per μm²
    nframes=2000,
    framerate=100.0,
    camera=camera
)

# Visualize with time coloring
fig = Figure(size=(800, 600))

ax = Axis(fig[1, 1], 
    title="Temporal Dynamics",
    xlabel="x (μm)",
    ylabel="y (μm)",
    aspect=DataAspect(),
    yreversed=true
)

# Get frame numbers and coordinates
frames = [e.frame for e in smld_noisy.emitters]
x = [e.x for e in smld_noisy.emitters]
y = [e.y for e in smld_noisy.emitters]

# Color by frame number
scatter!(ax, x, y, 
    color=frames, 
    colormap=:viridis, 
    markersize=4, 
    alpha=0.6
)

Colorbar(fig[1, 2], colormap=:viridis, label="Frame")

fig
```

## Creating SMLM Movies

This example shows how to create a frame-by-frame movie of SMLM data acquisition.

```julia
using SMLMSim
using CairoMakie

# Simulate data with circular patterns
pattern = Nmer2D(n=8, d=0.2)
camera = IdealCamera(1:256, 1:256, 0.1)
smld_true, smld_model, smld_noisy = simulate(
    pattern=pattern,
    ρ=1.0,
    nframes=500,
    framerate=50.0,
    camera=camera
)

# Create a function to visualize a specific frame
function visualize_frame(smld, frame_number)
    # Filter emitters in this frame
    frame_emitters = filter(e -> e.frame == frame_number, smld.emitters)
    
    # Extract coordinates
    x = [e.x for e in frame_emitters]
    y = [e.y for e in frame_emitters]
    photons = [e.photons for e in frame_emitters]
    
    # Create figure
    fig = Figure(size=(700, 600))
    
    ax = Axis(fig[1, 1], 
        title="Frame $frame_number",
        xlabel="x (μm)",
        ylabel="y (μm)",
        aspect=DataAspect(),
        yreversed=true
    )
    
    # Current frame emitters
    scatter!(ax, x, y, 
        color=photons,
        colormap=:plasma,
        markersize=8,
        alpha=0.8
    )
    
    # Also show accumulated emitters up to this frame
    acc_emitters = filter(e -> e.frame <= frame_number, smld.emitters)
    acc_x = [e.x for e in acc_emitters]
    acc_y = [e.y for e in acc_emitters]
    
    scatter!(ax, acc_x, acc_y, 
        color=:gray,
        markersize=3,
        alpha=0.2
    )
    
    Colorbar(fig[1, 2], colormap=:plasma, label="Photons")
    
    return fig
end

# Create movie (example for a specific frame)
frame_50 = visualize_frame(smld_noisy, 50)
```

For creating an actual movie, you would use the `CairoMakie.record` function:

```julia
# Code to create movie (not executed here)
frames = 1:10:500  # Subset of frames for the movie
record(fig, "smlm_movie.mp4", frames; framerate=10) do frame
    visualize_frame(smld_noisy, frame)
end
```