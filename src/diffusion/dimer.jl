
"""
    get_dimers(system::DiffusingMoleculeSystem)

Extract a new DiffusingMoleculeSystem containing only molecules in dimer state.

# Arguments
- `system::DiffusingMoleculeSystem`: Original system with all molecules

# Returns
- `DiffusingMoleculeSystem`: New system containing only dimers
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
"""
function get_dimers(systems::Vector{<:DiffusingMoleculeSystem})
    [get_dimers(system) for system in systems]
end

"""
    gen_dimer_images(systems::Vector{DiffusingMoleculeSystem}, psf::AbstractPSF; kwargs...)

Generate microscope images showing only dimers.

# Arguments
- `systems::Vector{DiffusingMoleculeSystem}`: Sequence of system states
- `psf::PSF`: Point spread function model

# Optional kwargs
- `frame_integration::Int=1`: Number of frames to integrate
- `poisson_noise::Bool=true`: Whether to add Poisson noise

# Returns
- `Array{Float64,3}`: Stack of dimer-only images [ny, nx, frames]
"""
function gen_dimer_images(systems::Vector{<:DiffusingMoleculeSystem}, psf::AbstractPSF;
                         frame_integration::Int=1, poisson_noise::Bool=true)
    dimer_systems = get_dimers(systems)
    gen_image_sequence(psf, dimer_systems; 
                      frame_integration=frame_integration,
                      poisson_noise=poisson_noise)
end

"""
    analyze_dimer_fraction(systems::Vector{DiffusingMoleculeSystem})

Calculate the fraction of molecules in dimer state over time.

# Arguments
- `systems::Vector{DiffusingMoleculeSystem}`: Sequence of system states

# Returns
- `Vector{Float64}`: Dimer fraction at each timepoint
"""
function analyze_dimer_fraction(systems::Vector{<:DiffusingMoleculeSystem})
    map(systems) do system
        n_dimers = count(mol -> mol.state == 2, system.molecules)
        n_total = length(system.molecules)
        n_dimers / n_total
    end
end

