"""
    Analysis functions for diffusion simulations.

This file contains functions for analyzing diffusion simulations,
including dimer detection and analysis of diffusion dynamics.
"""

"""
    get_dimers(system::DiffusingMoleculeSystem)

Extract a new DiffusingMoleculeSystem containing only molecules in dimer state.

# Arguments
- `system::DiffusingMoleculeSystem`: Original system with all molecules

# Returns
- `DiffusingMoleculeSystem`: New system containing only dimers

# Example
```julia
# Extract only dimers from a system
dimer_system = get_dimers(system)
```
"""
function get_dimers(system::DiffusingMoleculeSystem)
    # Extract molecules in dimer state
    dimer_molecules = filter(mol -> mol.state == 2, system.molecules)
    
    # Create new system with same parameters but only dimer molecules
    DiffusingMoleculeSystem(
        dimer_molecules,
        system.camera,
        system.box_size,
        system.n_frames,
        system.n_datasets,
        copy(system.metadata)
    )
end

"""
    get_dimers(systems::Vector{DiffusingMoleculeSystem})

Extract dimers from a sequence of system states.

# Arguments
- `systems::Vector{DiffusingMoleculeSystem}`: Sequence of system states

# Returns
- `Vector{DiffusingMoleculeSystem}`: Sequence of dimer-only states

# Example
```julia
# Run simulation
params = SmoluchowskiParams()
systems = simulate(params)

# Extract dimers at each timepoint
dimer_systems = get_dimers(systems)
```
"""
function get_dimers(systems::Vector{<:DiffusingMoleculeSystem})
    [get_dimers(system) for system in systems]
end



"""
    analyze_dimer_fraction(systems::Vector{DiffusingMoleculeSystem})

Calculate the fraction of molecules in dimer state over time.

# Arguments
- `systems::Vector{DiffusingMoleculeSystem}`: Sequence of system states

# Returns
- `Vector{Float64}`: Dimer fraction at each timepoint

# Example
```julia
# Run simulation
params = SmoluchowskiParams()
systems = simulate(params)

# Calculate dimer fraction over time
dimer_fractions = analyze_dimer_fraction(systems)

# Plot the results
using Plots
times = range(0, params.t_max, length=length(systems))
plot(times, dimer_fractions, 
     xlabel="Time (s)", 
     ylabel="Fraction of molecules in dimers",
     title="Dimer formation dynamics")
```
"""
function analyze_dimer_fraction(systems::Vector{<:DiffusingMoleculeSystem})
    map(systems) do system
        n_dimers = count(mol -> mol.state == 2, system.molecules)
        n_total = length(system.molecules)
        n_dimers / n_total
    end
end

"""
    analyze_dimer_lifetime(systems::Vector{DiffusingMoleculeSystem})

Analyze the lifetime distribution of dimers.

# Arguments
- `systems::Vector{DiffusingMoleculeSystem}`: Sequence of system states

# Returns
- `Vector{Float64}`: Distribution of dimer lifetimes

# Note
This function is not yet fully implemented.
"""
function analyze_dimer_lifetime(systems::Vector{<:DiffusingMoleculeSystem})
    # Get time step from metadata
    if !haskey(systems[1].metadata, "simulation_parameters")
        error("System metadata missing simulation parameters")
    end
    
    params = systems[1].metadata["simulation_parameters"]
    dt = params.dt
    
    # TODO: Implement dimer lifetime tracking
    @warn "Dimer lifetime analysis not yet fully implemented"
    
    return Float64[]
end