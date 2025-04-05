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



