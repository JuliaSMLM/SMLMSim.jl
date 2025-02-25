"""
    show_frame(system::DiffusingMoleculeSystem;
              title::String="", show_dimers::Bool=true)

Display a scatter plot of molecular positions in the system.

# Arguments
- `system::DiffusingMoleculeSystem`: System state to visualize

# Keyword Arguments
- `title::String=""`: Plot title
- `show_dimers::Bool=true`: Whether to color-code dimers differently

# Returns
- `Figure`: Generated figure

# Example
```julia
# Display the current state of a system
params = SmoluchowskiParams()
systems = simulate(params)
fig = show_frame(systems[50])  # Show the 50th time point
```
"""
function show_frame(system::DiffusingMoleculeSystem; 
                   title::String="",
                   show_dimers::Bool=true)
    
    # Extract positions
    x = [mol.x for mol in system.molecules]
    y = [mol.y for mol in system.molecules]
    
    # Determine colors
    colors = if show_dimers
        [mol.state == 2 ? :red : :blue for mol in system.molecules]
    else
        fill(:blue, length(system.molecules))
    end
    
    # Create figure
    f = Figure()
    ax = Axis(f[1,1]; 
        title = isempty(title) ? "Molecule Positions" : title,
        xlabel = "x (μm)",
        ylabel = "y (μm)",
        aspect = DataAspect(),
        yreversed = true)
    
    # Plot molecules
    scatter!(ax, x, y, 
        markersize = 10,
        color = colors)
    
    # Set limits based on box size
    xlims!(ax, (0, system.box_size))
    ylims!(ax, (system.box_size, 0))
    
    # Add legend if showing dimers
    if show_dimers
        n_monomers = count(mol -> mol.state == 1, system.molecules)
        n_dimers = count(mol -> mol.state == 2, system.molecules)
        
        Legend(f[1, 2],
            [MarkerElement(color=:blue, marker=:circle),
             MarkerElement(color=:red, marker=:circle)],
            ["Monomers ($n_monomers)", "Dimers ($n_dimers)"],
            "Molecule Types")
    end
    
    return f
end

"""
    show_frame(system::DiffusingMoleculeSystem, filename::String; kwargs...)

Display and save a scatter plot of molecular positions.

# Arguments
- `system::DiffusingMoleculeSystem`: System state to visualize
- `filename::String`: Output file path

# Keyword Arguments
- Passed to `show_frame(system; kwargs...)`

# Returns
- `Figure`: Generated figure

# Example
```julia
# Save the current state to a file
show_frame(systems[100], "snapshot_100.png"; 
          title="System state at t=1.0s")
```
"""
function show_frame(system::DiffusingMoleculeSystem, filename::String; kwargs...)
    f = show_frame(system; kwargs...)
    save(filename, f)
    return f
end

"""
    visualize_sequence(systems::Vector{DiffusingMoleculeSystem};
                      filename::String="smoluchowski.mp4",
                      framerate::Int=30,
                      show_dimers::Bool=true)

Generate an animation of molecular dynamics.

# Arguments
- `systems::Vector{DiffusingMoleculeSystem}`: Sequence of system states

# Keyword Arguments
- `filename::String="smoluchowski.mp4"`: Output file path
- `framerate::Int=30`: Video frame rate
- `show_dimers::Bool=true`: Whether to color-code dimers

# Returns
- `Nothing`: Saves animation to file

# Example
```julia
# Run simulation and create animation
params = SmoluchowskiParams(t_max=5.0)
systems = simulate(params)
visualize_sequence(systems, filename="diffusion_5s.mp4", framerate=60)
```
"""
function visualize_sequence(systems::Vector{<:DiffusingMoleculeSystem};
                          filename::String="smoluchowski.mp4",
                          framerate::Int=30,
                          show_dimers::Bool=true)
    
    # Input validation
    if isempty(systems)
        throw(ArgumentError("Systems vector cannot be empty"))
    end
    
    if framerate <= 0
        throw(ArgumentError("Frame rate must be positive"))
    end
    
    # Set up figure
    fig = Figure()
    ax = Axis(fig[1, 1], 
        xlabel = "x (μm)",
        ylabel = "y (μm)",
        title = "Molecular Dynamics Simulation",
        aspect = DataAspect(),
        yreversed = true)
    
    # Initial scatter plot
    x = [mol.x for mol in systems[1].molecules]
    y = [mol.y for mol in systems[1].molecules]
    colors = show_dimers ? 
        [mol.state == 2 ? :red : :blue for mol in systems[1].molecules] :
        fill(:blue, length(systems[1].molecules))
    
    sc = scatter!(ax, x, y, 
        color = colors,
        markersize = 10)
    
    # Set consistent limits
    box_size = systems[1].box_size
    xlims!(ax, (0, box_size))
    ylims!(ax, (box_size, 0))
    
    # Add legend if showing dimers
    if show_dimers
        leg = Legend(fig[1, 2], ["Monomers", "Dimers"], [:blue, :red])
    end
    
    # Information display
    info_text = Observable("Frame: 1")
    fig[2, :] = Label(fig, info_text)
    
    # Generate animation frames
    record(fig, filename, 1:length(systems); framerate=framerate) do i
        # Update title with frame info
        time = get(systems[i].metadata, "time", (i-1)*get(systems[i].metadata, "simulation_parameters", (dt=0.01,)).dt)
        ax.title = @sprintf("t = %.2f s", time)
        
        # Update positions and colors
        sc[1] = [mol.x for mol in systems[i].molecules]
        sc[2] = [mol.y for mol in systems[i].molecules]
        if show_dimers
            sc.color = [mol.state == 2 ? :red : :blue for mol in systems[i].molecules]
            
            # Update legend with counts
            n_monomers = count(mol -> mol.state == 1, systems[i].molecules)
            n_dimers = count(mol -> mol.state == 2, systems[i].molecules)
            leg.labels = ["Monomers ($n_monomers)", "Dimers ($n_dimers)"]
        end
        
        # Update info text
        info_text[] = @sprintf("Frame: %d/%d, Time: %.2f s", 
                               i, length(systems), time)
    end
    
    return nothing
end

"""
    visualize_simulation(params::SmoluchowskiParams;
                        filename::String="smoluchowski.mp4",
                        framerate::Int=30)

Run a simulation and generate visualization.

# Arguments
- `params::SmoluchowskiParams`: Simulation parameters

# Keyword Arguments
- `filename::String="smoluchowski.mp4"`: Output file path
- `framerate::Int=30`: Video frame rate

# Returns
- `Vector{DiffusingMoleculeSystem}`: Generated system states

# Example
```julia
# Set parameters and run simulation with visualization in one step
params = SmoluchowskiParams(
    density = 1.0,
    box_size = 20.0,
    t_max = 10.0
)
systems = visualize_simulation(params, filename="high_density.mp4")
```
"""
function visualize_simulation(params::SmoluchowskiParams;
                            filename::String="smoluchowski.mp4",
                            framerate::Int=30)
    
    # Run simulation
    systems = simulate(params)
    
    # Generate visualization
    visualize_sequence(systems; 
        filename=filename,
        framerate=framerate)
    
    return systems
end