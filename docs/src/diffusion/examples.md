```@meta
CurrentModule = SMLMSim
```

# Diffusion-Interaction Examples

This page provides complete examples for using the diffusion-interaction simulation capabilities of SMLMSim.

## Working with Trajectories

SMLMSim provides several utility functions for working with trajectories:

- `get_track(smld, id)`: Returns a new SMLD containing only emitters with the specified track_id
- `get_num_tracks(smld)`: Returns the number of unique tracks in an SMLD
- `get_tracks(smld)`: Returns a vector of SMDLs, one for each unique track

These functions are useful for analyzing and visualizing trajectory data, as demonstrated in the examples below.

## Basic Diffusion Simulation

This example demonstrates how to run a basic diffusion simulation and visualize the results:

```@example
using SMLMSim
using CairoMakie
using MicroscopePSFs

# Set up simulation parameters
params = DiffusionSMLMParams(
    density = 2.0,        # molecules per μm²
    box_size = 10.0,      # μm
    diff_monomer = 0.1,   # μm²/s
    diff_dimer = 0.05,    # μm²/s
    k_off = 0.2,          # s⁻¹
    r_react = 0.01,       # μm
    d_dimer = 0.05,       # μm
    dt = 0.01,            # s
    t_max = 10.0,         # s
    camera_framerate = 10.0  # fps
)

# Run simulation - returns (smld, SimInfo) tuple
smld, sim_info = simulate(params; photons=1000.0)

# Extract coordinates based on monomer/dimer state for a specific frame
function extract_frame_by_state(smld, frame_num)
    # This approach doesn't use our new utility functions directly yet,
    # but demonstrates what would be useful additions in the future
    # (e.g., get_frame and filter_by_state functions)
    
    # Filter emitters for this frame the traditional way
    frame_emitters = filter(e -> e.frame == frame_num, smld.emitters)
    
    # Extract coordinates and states
    x = [e.x for e in frame_emitters]
    y = [e.y for e in frame_emitters]
    states = [e.state for e in frame_emitters]
    
    # Create state-based colors
    colors = map(states) do state
        state == :monomer ? :blue : :red
    end
    
    return x, y, colors
end

# Create visualization for multiple frames
frames_to_show = [1, 20, 40, 60, 80, 100]
fig = Figure(size=(1000, 600))

for (i, frame) in enumerate(frames_to_show)
    # Calculate subplot position (2x3 grid)
    row = div(i-1, 3) + 1
    col = rem(i-1, 3) + 1
    
    # Create axis
    ax = Axis(fig[row, col], 
        title="Frame $frame",
        xlabel="x (μm)",
        ylabel="y (μm)",
        aspect=DataAspect()
    )
    
    # Extract and plot data
    x, y, colors = extract_frame_by_state(smld, frame)
    scatter!(ax, x, y, color=colors, markersize=6)
    
    # Set consistent axis limits
    limits!(ax, 0, params.box_size, 0, params.box_size)
end

# Add legend
fig[1:2, 4] = Legend(fig, 
    [MarkerElement(color=:blue, marker=:circle), 
     MarkerElement(color=:red, marker=:circle)],
    ["Monomer", "Dimer"],
    "Molecular States"
)

fig

```

## Frame Integration for Time-Lapse Imaging

This example demonstrates how to use frame integration to create realistic time-lapse microscopy data from diffusion simulations:

```@example
using SMLMSim
using MicroscopePSFs
using CairoMakie

# Set up diffusion parameters with high temporal resolution
params = DiffusionSMLMParams(
    density = 1.0,          # molecules per μm²
    box_size = 10.0,        # μm
    diff_monomer = 0.2,     # μm²/s (moderate diffusion)
    diff_dimer = 0.1,       # μm²/s
    k_off = 0.1,            # s⁻¹
    r_react = 0.01,         # μm
    d_dimer = 0.05,         # μm
    dt = 0.001,             # s (1ms time steps for accurate diffusion)
    t_max = 5.0,            # s (5 seconds total)
    camera_framerate = 10.0, # fps (10 frames per second)
    camera_exposure = 0.1    # s (100ms exposure time)
)

# Run simulation - returns (smld, SimInfo) tuple
smld, sim_info = simulate(params)

# Setup camera and PSF
pixelsize = 0.1  # 100nm pixels
pixels = Int64(round(params.box_size/pixelsize))
camera = IdealCamera(1:pixels, 1:pixels, pixelsize)

# Create PSF model
psf = MicroscopePSFs.GaussianPSF(0.15)  # 150nm PSF width

# Generate microscope images from simulation - returns (images, ImageInfo) tuple
# For diffusion simulations, the camera integration time (exposure) has already been
# modeled in the simulation process, so each frame already includes the positions
# from all emitters that appeared during the exposure window
images, img_info = gen_images(smld, psf;
    bg=5.0,               # background photons per pixel
    poisson_noise=true     # add photon counting noise
)

# Display a single frame
fig = Figure(size=(600, 600))
ax = Axis(fig[1, 1],
    title="Diffusion with Frame Integration",
    aspect=DataAspect(),
    yreversed=true
)

# Select a frame to display
frame_to_show = 3  # Changed from 15 to a valid frame index (between 1 and 6)
heatmap!(ax, transpose(images[:, :, frame_to_show]), colormap=:inferno)

fig
```

The simulation already handles motion blur effects in a realistic way:

- The `camera_exposure` parameter in the simulation determines how long each camera frame integrates photons
- During the exposure window, multiple emitter positions from the same track_id are captured
- This naturally creates motion blur effects where fast-moving particles appear more blurred
- The resulting images accurately represent what would be seen in real microscopy experiments

## Analyzing Dimer Formation

This example demonstrates how to analyze dimer formation dynamics:

```@example
using SMLMSim
using CairoMakie

# Set up simulations with different dissociation rates
params_stable = DiffusionSMLMParams(
    density = 1.0,        # molecules per μm²
    diff_monomer = 0.1,   # μm²/s
    diff_dimer = 0.05,    # μm²/s
    k_off = 0.05,         # s⁻¹ (stable dimers)
    t_max = 20.0          # s
)

params_unstable = DiffusionSMLMParams(
    density = 1.0,        # molecules per μm²
    diff_monomer = 0.1,   # μm²/s
    diff_dimer = 0.05,    # μm²/s
    k_off = 0.5,          # s⁻¹ (unstable dimers)
    t_max = 20.0          # s
)

# Run simulations - returns (smld, SimInfo) tuples
smld_stable, _ = simulate(params_stable)
smld_unstable, _ = simulate(params_unstable)

# Analyze dimer formation
frames_stable, frac_stable = analyze_dimer_fraction(smld_stable)
frames_unstable, frac_unstable = analyze_dimer_fraction(smld_unstable)

# Calculate time in seconds
time_stable = (frames_stable .- 1) ./ params_stable.camera_framerate
time_unstable = (frames_unstable .- 1) ./ params_unstable.camera_framerate

# Visualize dimer formation dynamics
fig = Figure(size=(800, 500))

ax = Axis(fig[1, 1],
    title="Dimer Formation Dynamics",
    xlabel="Time (s)",
    ylabel="Fraction of molecules in dimers"
)

lines!(ax, time_stable, frac_stable, linewidth=3, color=:blue, 
       label="Stable (k_off=0.05 s⁻¹)")
lines!(ax, time_unstable, frac_unstable, linewidth=3, color=:red, 
       label="Unstable (k_off=0.5 s⁻¹)")

axislegend(ax)

fig

```

## Generating Microscope Images

This example shows how to generate microscope images from diffusion simulations:

```@example
using SMLMSim
using CairoMakie
using MicroscopePSFs

# Set simulation parameters
params = DiffusionSMLMParams(
    density = 0.3,        # molecules per μm²
    box_size = 10.0,      # μm
    diff_monomer = 0.1,   # μm²/s
    diff_dimer = 0.05,    # μm²/s
    k_off = 0.2,          # s⁻¹
    r_react = 0.01,       # μm
    d_dimer = 0.05,       # μm
    dt = 0.01,            # s
    t_max = 5.0,          # s
    camera_framerate = 10.0,  # fps
    camera_exposure = 0.1     # s
)

# Run simulation - returns (smld, SimInfo) tuple
smld, sim_info = simulate(params)

# Set up camera and PSF
pixelsize = 0.1  # 100nm pixels
pixels = Int64(round(params.box_size/pixelsize))
camera = IdealCamera(1:pixels, 1:pixels, pixelsize)

# Set up PSF (Gaussian with 150nm width)
psf = MicroscopePSFs.GaussianPSF(0.15)  # 150nm PSF width

# Generate images for all molecules - returns (images, ImageInfo) tuple
images_all, _ = gen_images(smld, psf;
    bg=5.0,
    poisson_noise=true
)

# Extract only dimers
dimer_smld = get_dimers(smld)

# Generate images showing only dimers - returns (images, ImageInfo) tuple
images_dimers, _ = gen_images(dimer_smld, psf;
    bg=5.0,
    poisson_noise=true
)

# Visualize images
function display_frames(images_all, images_dimers, frame_indices)
    fig = Figure(size=(1200, 800))
    
    for (i, frame) in enumerate(frame_indices)
        # Row 1: All molecules
        ax1 = Axis(fig[1, i], 
            title="All Molecules - Frame $frame",
            aspect=DataAspect(),
            yreversed=true
        )
        hidedecorations!(ax1)
        
        # Row 2: Dimers only
        ax2 = Axis(fig[2, i], 
            title="Dimers Only - Frame $frame",
            aspect=DataAspect(),
            yreversed=true
        )
        hidedecorations!(ax2)
        
        # Display images
        heatmap!(ax1, transpose(images_all[:, :, frame]), colormap=:inferno)
        heatmap!(ax2, transpose(images_dimers[:, :, frame]), colormap=:inferno)
    end
    
    return fig
end

# Show frames 10, 20, 30
frame_indices = [10, 20, 30]
fig = display_frames(images_all, images_dimers, frame_indices)
```

## Two Interacting Particles

This example shows two particles interacting in a small box with reflecting boundary conditions, using the new starting conditions feature to place them at specific initial positions:

```@example
using SMLMSim
using CairoMakie

# Set up a minimal simulation with just two particles
params = DiffusionSMLMParams(
    box_size = 1.0,       # 1 μm box for close interactions
    diff_monomer = 0.1,   # μm²/s
    diff_dimer = 0.05,    # μm²/s
    k_off = 0.5,          # s⁻¹ (moderate dimer stability)
    r_react = 0.05,       # μm (large reaction radius for demonstration)
    d_dimer = 0.07,       # μm (dimer separation)
    dt = 0.01,            # s
    t_max = 5.0,          # s
    boundary = "reflecting",  # Reflecting boundaries
    camera_framerate = 10.0   # fps
)

# Create two particles with specific initial positions
particle1 = DiffusingEmitter2D{Float64}(
    0.2, 0.2,       # Position in lower-left quadrant
    1000.0,         # Photons
    0.0,            # Initial timestamp
    1,              # Initial frame
    1,              # Dataset
    1,              # track_id
    :monomer,       # Initial state
    nothing         # No partner initially
)

particle2 = DiffusingEmitter2D{Float64}(
    0.8, 0.8,       # Position in upper-right quadrant
    1000.0,         # Photons
    0.0,            # Initial timestamp
    1,              # Initial frame
    1,              # Dataset
    2,              # track_id
    :monomer,       # Initial state
    nothing         # No partner initially
)

# Run simulation with custom starting positions - returns (smld, SimInfo) tuple
smld, sim_info = simulate(params; starting_conditions=[particle1, particle2])

track_smlds = get_tracks(smld)

# Convert to the format needed for plotting
trajectories = []
for track_smld in track_smlds
    # Get ID from first emitter
    id = track_smld.emitters[1].track_id
    
    # Sort by timestamp
    sort!(track_smld.emitters, by = e -> e.timestamp)
    
    # Extract coordinates and state
    times = [e.timestamp for e in track_smld.emitters]
    x = [e.x for e in track_smld.emitters]
    y = [e.y for e in track_smld.emitters]
    states = [e.state for e in track_smld.emitters]
    
    push!(trajectories, (id=id, times=times, x=x, y=y, states=states))
end

# Visualize interaction dynamics
fig = Figure(size=(700, 600))

ax = Axis(fig[1, 1], 
    title="Two Particles in 1μm Box (Reflecting Boundaries)",
    xlabel="x (μm)",
    ylabel="y (μm)",
    aspect=DataAspect()
)

# Plot trajectories with state-dependent coloring
for (i, traj) in enumerate(trajectories)
    # Create segments with colors based on state
    segments_x = []
    segments_y = []
    colors = []
    
    for j in 1:(length(traj.times)-1)
        push!(segments_x, [traj.x[j], traj.x[j+1]])
        push!(segments_y, [traj.y[j], traj.y[j+1]])
        push!(colors, traj.states[j] == :monomer ? :blue : :red)
    end
    
    # Plot each segment with appropriate color
    for j in 1:length(segments_x)
        lines!(ax, segments_x[j], segments_y[j], 
               color=colors[j], linewidth=2, 
               label=j==1 ? "Particle $(traj.id)" : nothing)
    end
    
    # Mark starting position
    scatter!(ax, [traj.x[1]], [traj.y[1]], 
            color=:black, marker=:circle, markersize=10)
    
    # Mark ending position
    scatter!(ax, [traj.x[end]], [traj.y[end]], 
            color=:black, marker=:star, markersize=12)
end

# Show box boundaries
box = [0 0; 1 0; 1 1; 0 1; 0 0]
lines!(ax, box[:, 1], box[:, 2], color=:black, linewidth=2)

# Add legend for state colors
legend_elements = [
    LineElement(color=:blue, linewidth=3),
    LineElement(color=:red, linewidth=3),
    MarkerElement(color=:black, marker=:circle, markersize=8),
    MarkerElement(color=:black, marker=:star, markersize=10)
]
legend_labels = ["Monomer", "Dimer", "Start", "End"]

Legend(fig[1, 2], legend_elements, legend_labels, "States")

# Set axis limits with some padding
limits!(ax, -0.05, 1.05, -0.05, 1.05)

fig

```

