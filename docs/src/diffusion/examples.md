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
    density = 0.5,        # molecules per μm²
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

# Run simulation
smld = simulate(params; photons=1000.0)

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

# Run simulation
smld = simulate(params)

# Setup camera and PSF
pixelsize = 0.1  # 100nm pixels
pixels = Int64(round(params.box_size/pixelsize))
camera = IdealCamera(1:pixels, 1:pixels, pixelsize)

# Create PSF model
psf = MicroscopePSFs.GaussianPSF(0.15)  # 150nm PSF width

# Generate microscope images with frame integration
# The frame_integration parameter determines how many simulation time
# points are integrated into each output frame
# Note: For diffusion simulations, the smld already contains correct positions without 
# localization uncertainty, so we can use it directly for generating camera images
images = gen_images(smld, psf;
    bg=5.0,               # background photons per pixel
    frame_integration=10,  # integrate 10 simulation steps per frame
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

The `frame_integration` parameter is crucial for realistic diffusion imaging:

- With high values, each camera frame integrates multiple diffusion steps, creating motion blur effects typical in real microscopy
- The integrated frames combine positions from multiple simulation time points, accurately modeling camera exposure
- Fast-moving particles appear more blurred while slow or stationary ones remain sharp

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

# Run simulations
smld_stable = simulate(params_stable)
smld_unstable = simulate(params_unstable)

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

# Run simulation
smld = simulate(params)

# Set up camera and PSF
pixelsize = 0.1  # 100nm pixels
pixels = Int64(round(params.box_size/pixelsize))
camera = IdealCamera(1:pixels, 1:pixels, pixelsize)

# Set up PSF (Gaussian with 150nm width)
psf = MicroscopePSFs.GaussianPSF(0.15)  # 150nm PSF width

# Generate images for all molecules
images_all = gen_images(smld, psf;
    bg=5.0,
    poisson_noise=true
)

# Extract only dimers
dimer_smld = get_dimers(smld)

# Generate images showing only dimers
images_dimers = gen_images(dimer_smld, psf;
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

# Run simulation with custom starting positions
smld = simulate(params; starting_conditions=[particle1, particle2])

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

## Multi-Stage Simulation with Modified Parameters

This example demonstrates how to run a simulation in multiple stages with different parameters, using the final state of one simulation as the starting point for the next:

```@example
using SMLMSim
using CairoMakie

# Stage 1: Initial diffusion with specific starting positions
params1 = DiffusionSMLMParams(
    box_size = 1.0,       # 1 μm box
    diff_monomer = 0.05,  # μm²/s (slow diffusion)
    diff_dimer = 0.02,    # μm²/s
    k_off = 0.2,          # s⁻¹
    r_react = 0.05,       # μm
    d_dimer = 0.07,       # μm
    dt = 0.01,            # s
    t_max = 2.0,          # s (short first stage)
    boundary = "reflecting",
    camera_framerate = 10.0
)

# Create initial positions for particles at opposite corners
particle1 = DiffusingEmitter2D{Float64}(0.1, 0.1, 1000.0, 0.0, 1, 1, 1, :monomer, nothing)
particle2 = DiffusingEmitter2D{Float64}(0.9, 0.9, 1000.0, 0.0, 1, 1, 2, :monomer, nothing)

# Run first stage of simulation
smld1 = simulate(params1; starting_conditions=[particle1, particle2])

# Extract the final state from the first simulation
final_state = extract_final_state(smld1)

# Stage 2: Continue with faster diffusion and different dissociation rate
params2 = DiffusionSMLMParams(
    box_size = 1.0,       # Same box size
    diff_monomer = 0.2,   # μm²/s (faster diffusion)
    diff_dimer = 0.1,     # μm²/s
    k_off = 0.05,         # s⁻¹ (more stable dimers)
    r_react = 0.05,       # μm
    d_dimer = 0.07,       # μm
    dt = 0.01,            # s
    t_max = 3.0,          # s (additional 3 seconds)
    boundary = "reflecting",
    camera_framerate = 10.0
)

# Run second stage using final state from first stage
smld2 = simulate(params2; starting_conditions=final_state)

# Combine results from both stages to show full trajectory
# Adjust timestamps and frames for stage 2 to continue from stage 1
combined_emitters = []

# Add all emitters from stage 1
append!(combined_emitters, smld1.emitters)

# Add emitters from stage 2 with adjusted timestamps and frames
max_time = maximum([e.timestamp for e in smld1.emitters])
max_frame = maximum([e.frame for e in smld1.emitters])

for e in smld2.emitters
    if isa(e, DiffusingEmitter2D)
        # Create new emitter with adjusted timestamp and frame
        adjusted_emitter = DiffusingEmitter2D{typeof(e.x)}(
            e.x, e.y,                 # Position
            e.photons,                # Photons
            e.timestamp + max_time,   # Adjusted timestamp
            e.frame + max_frame,      # Adjusted frame
            e.dataset,                # Dataset
            e.id,                     # ID
            e.state,                  # State
            e.partner_id              # Partner ID
        )
        push!(combined_emitters, adjusted_emitter)
    end
end

# Create combined SMLD
combined_metadata = copy(smld1.metadata)
combined_metadata["simulation_type"] = "multi_stage"
combined_metadata["stage1_params"] = params1
combined_metadata["stage2_params"] = params2

camera = IdealCamera(1:10, 1:10, 0.1)  # Simple camera for visualization
combined_smld = BasicSMLD(
    combined_emitters,
    camera,
    max_frame + maximum([e.frame for e in smld2.emitters]),
    1,
    combined_metadata
)

# Get trajectories from combined data
track_smlds = get_tracks(combined_smld)

# Process trajectories for visualization
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

# Create visualization showing both stages
fig = Figure(size=(900, 400))

# Create two side-by-side plots
ax1 = Axis(fig[1, 1], 
    title="Stage 1: Slow Diffusion",
    xlabel="x (μm)",
    ylabel="y (μm)",
    aspect=DataAspect()
)

ax2 = Axis(fig[1, 2], 
    title="Stage 2: Fast Diffusion",
    xlabel="x (μm)",
    ylabel="y (μm)",
    aspect=DataAspect()
)

# Box boundaries
box = [0 0; 1 0; 1 1; 0 1; 0 0]
lines!(ax1, box[:, 1], box[:, 2], color=:black, linewidth=1)
lines!(ax2, box[:, 1], box[:, 2], color=:black, linewidth=1)

# Plot trajectories with different colors for each stage
for (i, traj) in enumerate(trajectories)
    # Determine the time cutoff between stages
    stage_separation = max_time
    
    # Find indices for each stage
    stage1_indices = findall(t -> t <= stage_separation, traj.times)
    stage2_indices = findall(t -> t > stage_separation, traj.times)
    
    # Plot stage 1 trajectory
    if !isempty(stage1_indices)
        x1 = traj.x[stage1_indices]
        y1 = traj.y[stage1_indices]
        states1 = traj.states[stage1_indices]
        
        # Plot segments with state-dependent colors
        for j in 1:(length(stage1_indices)-1)
            lines!(ax1, 
                [x1[j], x1[j+1]], 
                [y1[j], y1[j+1]], 
                color = states1[j] == :monomer ? :blue : :red,
                linewidth = 2
            )
        end
        
        # Mark start and end points
        scatter!(ax1, [x1[1]], [y1[1]], color=:black, marker=:circle, markersize=8)
        scatter!(ax1, [x1[end]], [y1[end]], color=:black, marker=:star, markersize=10)
    end
    
    # Plot stage 2 trajectory
    if !isempty(stage2_indices)
        x2 = traj.x[stage2_indices]
        y2 = traj.y[stage2_indices]
        states2 = traj.states[stage2_indices]
        
        # Plot segments with state-dependent colors
        for j in 1:(length(stage2_indices)-1)
            lines!(ax2, 
                [x2[j], x2[j+1]], 
                [y2[j], y2[j+1]], 
                color = states2[j] == :monomer ? :blue : :red,
                linewidth = 2
            )
        end
        
        # Mark start and end points
        scatter!(ax2, [x2[1]], [y2[1]], color=:black, marker=:circle, markersize=8)
        scatter!(ax2, [x2[end]], [y2[end]], color=:black, marker=:star, markersize=10)
    end
end

# Add a legend
fig[1, 3] = Legend(fig,
    [
        LineElement(color=:blue, linewidth=3),
        LineElement(color=:red, linewidth=3),
        MarkerElement(color=:black, marker=:circle, markersize=8),
        MarkerElement(color=:black, marker=:star, markersize=10)
    ],
    ["Monomer", "Dimer", "Start", "End"],
    "States"
)

# Set axis limits with padding
limits!(ax1, -0.05, 1.05, -0.05, 1.05)
limits!(ax2, -0.05, 1.05, -0.05, 1.05)

fig
```

This example demonstrates how the starting conditions feature allows for flexible multi-stage simulations where parameters like diffusion coefficients and reaction rates can be changed between stages. This can be useful for modeling complex scenarios such as:

- Environmental changes affecting molecular behavior
- Dynamic transitions between different experimental conditions
- Studying the effects of sudden parameter changes on molecular interactions
- Creating more complex and realistic biological models
