```@meta
CurrentModule = SMLMSim
```

# Diffusion-Interaction Examples

This page provides complete examples for using the diffusion-interaction simulation capabilities of SMLMSim.

## Basic Diffusion Simulation

This example demonstrates how to run a basic diffusion simulation and visualize the results:

```julia
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
smld = simulate(params)

# Extract coordinates from different frames
function extract_frame(smld, frame_num)
    # Filter emitters for this frame
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
    x, y, colors = extract_frame(smld, frame)
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

```julia
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
psf = MicroscopePSFs.Gaussian2D(0.15)  # 150nm PSF width

# Generate microscope images with frame integration
# The frame_integration parameter determines how many simulation time
# points are integrated into each output frame
images = gen_images(smld, psf;
    photons=2000.0,       # photons per emitter
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
frame_to_show = 15
heatmap!(ax, transpose(images[:, :, frame_to_show]), colormap=:inferno)

fig
```

The `frame_integration` parameter is crucial for realistic diffusion imaging:

- With high values, each camera frame integrates multiple diffusion steps, creating motion blur effects typical in real microscopy
- The integrated frames combine positions from multiple simulation time points, accurately modeling camera exposure
- Fast-moving particles appear more blurred while slow or stationary ones remain sharp

## Analyzing Dimer Formation

This example demonstrates how to analyze dimer formation dynamics:

```julia
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

axislegend(ax, position=:right)

fig
```

## Generating Microscope Images

This example shows how to generate microscope images from diffusion simulations:

```julia
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
psf = MicroscopePSFs.Gaussian2D(0.15)  # 150nm PSF width

# Generate images for all molecules
images_all = gen_images(smld, psf;
    photons=1000.0,
    bg=5.0,
    poisson_noise=true
)

# Extract only dimers
dimer_smld = get_dimers(smld)

# Generate images showing only dimers
images_dimers = gen_images(dimer_smld, psf;
    photons=1000.0,
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
    }
    
    return fig
end

# Show frames 10, 20, 30
frame_indices = [10, 20, 30]
display_frames(images_all, images_dimers, frame_indices)
```

## Diffusion with Different Boundary Conditions

This example compares periodic and reflecting boundary conditions:

```julia
using SMLMSim
using CairoMakie

# Set up simulations with different boundary conditions
params_periodic = DiffusionSMLMParams(
    density = 0.3,        # molecules per μm²
    box_size = 10.0,      # μm
    diff_monomer = 0.2,   # μm²/s (faster diffusion)
    t_max = 10.0,         # s
    boundary = "periodic"  # Periodic boundaries
)

params_reflecting = DiffusionSMLMParams(
    density = 0.3,        # molecules per μm²
    box_size = 10.0,      # μm
    diff_monomer = 0.2,   # μm²/s (faster diffusion)
    t_max = 10.0,         # s
    boundary = "reflecting"  # Reflecting boundaries
)

# Run simulations
smld_periodic = simulate(params_periodic)
smld_reflecting = simulate(params_reflecting)

# Function to extract trajectories
function extract_trajectories(smld)
    # Group emitters by ID
    emitters_by_id = Dict()
    for e in smld.emitters
        if !haskey(emitters_by_id, e.id)
            emitters_by_id[e.id] = []
        end
        push!(emitters_by_id[e.id], e)
    end
    
    # Convert to trajectories
    trajectories = []
    for (id, emitters) in emitters_by_id
        # Sort by timestamp
        sort!(emitters, by = e -> e.timestamp)
        
        # Extract coordinates
        times = [e.timestamp for e in emitters]
        x = [e.x for e in emitters]
        y = [e.y for e in emitters]
        
        push!(trajectories, (id=id, times=times, x=x, y=y))
    end
    
    return trajectories
end

# Get trajectories
traj_periodic = extract_trajectories(smld_periodic)
traj_reflecting = extract_trajectories(smld_reflecting)

# Visualize trajectories
function plot_trajectories(trajectories, title, box_size)
    fig = Figure(size=(700, 600))
    
    ax = Axis(fig[1, 1], 
        title=title,
        xlabel="x (μm)",
        ylabel="y (μm)",
        aspect=DataAspect()
    )
    
    # Plot first 20 trajectories
    for i in 1:min(20, length(trajectories))
        traj = trajectories[i]
        lines!(ax, traj.x, traj.y, color=Cycled(i), linewidth=1.5, alpha=0.7)
        
        # Mark start position
        scatter!(ax, [traj.x[1]], [traj.y[1]], color=Cycled(i), marker=:circle, 
                markersize=8)
        
        # Mark end position
        scatter!(ax, [traj.x[end]], [traj.y[end]], color=Cycled(i), marker=:star, 
                markersize=10)
    end
    
    # Show box boundaries
    box = [0 0; box_size 0; box_size box_size; 0 box_size; 0 0]
    lines!(ax, box[:, 1], box[:, 2], color=:black, linewidth=2)
    
    # Set axis limits with some padding
    limits!(ax, -0.5, box_size+0.5, -0.5, box_size+0.5)
    
    return fig
end

# Plot both types of trajectories
fig_periodic = plot_trajectories(traj_periodic, "Periodic Boundaries", 
                                params_periodic.box_size)
fig_reflecting = plot_trajectories(traj_reflecting, "Reflecting Boundaries", 
                                  params_reflecting.box_size)

# Display figures
fig_periodic, fig_reflecting
```

## Long-term Evolution of Dimer Population

This example simulates the long-term evolution of dimer formation under different conditions:

```julia
using SMLMSim
using CairoMakie

# Set up parameters for longer simulation
params = DiffusionSMLMParams(
    density = 0.5,        # molecules per μm²
    box_size = 20.0,      # μm
    diff_monomer = 0.1,   # μm²/s
    diff_dimer = 0.05,    # μm²/s
    k_off = 0.1,          # s⁻¹
    r_react = 0.01,       # μm
    d_dimer = 0.05,       # μm
    dt = 0.01,            # s
    t_max = 60.0,         # s (longer simulation)
    camera_framerate = 5.0 # fps (slower framerate for long simulation)
)

# Run simulation
smld = simulate(params)

# Calculate dimer fraction over time
frames, dimer_fractions = analyze_dimer_fraction(smld)

# Convert frames to time
time = (frames .- 1) ./ params.camera_framerate

# Calculate theoretical equilibrium dimer fraction
# For simple second-order kinetics with monomer-monomer association
# and first-order dimer dissociation
function theoretical_dimer_fraction(t, k_on, k_off, initial_concentration)
    # k_on is in μm²/s units
    # k_off is in s⁻¹ units
    # initial_concentration is in molecules/μm² units
    
    # Equilibrium constant (dimensionless)
    K_eq = k_on / k_off
    
    # Total concentration (constant)
    c_total = initial_concentration
    
    # Equilibrium dimer fraction
    f_eq = (1 + 4*K_eq*c_total - sqrt(1 + 8*K_eq*c_total)) / (4*K_eq*c_total)
    
    # Time-dependent approach to equilibrium
    # This is a simplified model assuming well-mixed conditions
    τ = 1 / (k_off + k_on * c_total)
    f_t = f_eq * (1 - exp(-t/τ))
    
    return f_t
end

# Estimate k_on from r_react and diffusion
k_on = 2π * params.diff_monomer * params.r_react
k_off = params.k_off
c_total = params.density

# Calculate theoretical curve
t_theory = LinRange(0, maximum(time), 1000)
f_theory = theoretical_dimer_fraction.(t_theory, k_on, k_off, c_total)

# Visualize results
fig = Figure(size=(800, 500))

ax = Axis(fig[1, 1],
    title="Dimer Formation Dynamics",
    xlabel="Time (s)",
    ylabel="Fraction of molecules in dimers"
)

# Plot simulation results
scatter!(ax, time, dimer_fractions, color=:blue, markersize=6, 
        label="Simulation")

# Plot theoretical curve
lines!(ax, t_theory, f_theory, color=:red, linewidth=2, 
      label="Theoretical model")

# Add equilibrium line
hlines!(ax, theoretical_dimer_fraction(Inf, k_on, k_off, c_total), 
       color=:black, linestyle=:dash, linewidth=1, 
       label="Equilibrium")

axislegend(ax, position=:right)

fig
```