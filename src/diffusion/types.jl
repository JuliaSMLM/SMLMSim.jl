"""
    DiffusingMolecule{E<:AbstractEmitter}

A molecule that can diffuse and form dimers, containing an SMLM emitter.

# Fields
- `emitter::E`: underlying SMLM emitter (2D/3D, fitted/unfitted)
- `state::Int`: oligomerization state (1=monomer, 2=dimer)
- `id::Int`: unique identifier for tracking
- `link::Union{Nothing,Int}`: id of linked molecule (for dimers)
- `updated::Bool`: flag for position updates

# Type Parameters
- `E`: Type of emitter contained in the molecule (must be a subtype of AbstractEmitter)

# Examples
```julia
# Create a 2D diffusing molecule
emitter = Emitter2D(1.0, 2.0, 1000.0)
molecule = DiffusingMolecule(emitter, 1, 1, nothing, false)

# Create a 3D diffusing molecule
emitter3d = Emitter3D(1.0, 2.0, 0.5, 1000.0)
molecule3d = DiffusingMolecule(emitter3d, 1, 2, nothing, false)
```
"""
mutable struct DiffusingMolecule{E<:AbstractEmitter}
    emitter::E
    state::Int  # 1=monomer, 2=dimer
    id::Int     # unique identifier
    link::Union{Nothing,Int}  # id of linked molecule
    updated::Bool
end

# Property access - handle both emitter fields and direct fields
function Base.getproperty(m::DiffusingMolecule, s::Symbol)
    if s ∈ (:emitter, :state, :id, :link, :updated)
        return getfield(m, s)
    else
        return getfield(m.emitter, s)
    end
end

function Base.setproperty!(m::DiffusingMolecule, s::Symbol, v)
    if s ∈ (:emitter, :state, :id, :link, :updated)
        setfield!(m, s, v)
    else
        setfield!(m.emitter, s, v)
    end
end

"""
    DiffusingMoleculeSystem{E<:AbstractEmitter} <: SMLD

System of diffusing molecules built on SMLD.

# Fields
- `molecules::Vector{DiffusingMolecule{E}}`: Vector of diffusing molecules
- `camera::AbstractCamera`: Camera used for acquisition
- `box_size::Float64`: Size of simulation box in microns
- `n_frames::Int`: Total number of frames
- `n_datasets::Int`: Number of datasets (typically 1)
- `metadata::Dict{String,Any}`: Additional simulation parameters

# Type Parameters
- `E`: Type of emitter contained in the molecules (must be a subtype of AbstractEmitter)

# Examples
```julia
# Create a simple system with a single molecule
camera = IdealCamera(1:100, 1:100, 0.1)
emitter = Emitter2D(1.0, 2.0, 1000.0)
molecule = DiffusingMolecule(emitter, 1, 1, nothing, false)
system = DiffusingMoleculeSystem(
    [molecule], 
    camera, 
    10.0,  # box size
    1,     # frames
    1,     # datasets
    Dict{String,Any}("simulation_type" => "diffusion")
)
```
"""
mutable struct DiffusingMoleculeSystem{E<:AbstractEmitter} <: SMLD
    molecules::Vector{DiffusingMolecule{E}}
    camera::SMLMData.AbstractCamera
    box_size::Float64
    n_frames::Int
    n_datasets::Int
    metadata::Dict{String,Any}
end

# Property access for DiffusingMoleculeSystem
function Base.getproperty(sys::DiffusingMoleculeSystem, s::Symbol)
    if s === :emitters
        return [m.emitter for m in getfield(sys, :molecules)]
    else
        return getfield(sys, s)
    end
end

# Helper functions for geometric calculations
"""
    calc_r(mol1::DiffusingMolecule{<:Emitter2D}, mol2::DiffusingMolecule{<:Emitter2D})

Calculate Euclidean distance between two 2D molecules.

# Arguments
- `mol1::DiffusingMolecule{<:Emitter2D}`: First molecule
- `mol2::DiffusingMolecule{<:Emitter2D}`: Second molecule

# Returns
- `Float64`: Distance between molecules in microns

# Example
```julia
distance = calc_r(molecule1, molecule2)
```
"""
function calc_r(mol1::DiffusingMolecule{<:Emitter2D}, mol2::DiffusingMolecule{<:Emitter2D})
    sqrt((mol1.x - mol2.x)^2 + (mol1.y - mol2.y)^2)
end

"""
    calc_r(mol1::DiffusingMolecule{<:Emitter3D}, mol2::DiffusingMolecule{<:Emitter3D})

Calculate Euclidean distance between two 3D molecules.

# Arguments
- `mol1::DiffusingMolecule{<:Emitter3D}`: First molecule
- `mol2::DiffusingMolecule{<:Emitter3D}`: Second molecule

# Returns
- `Float64`: Distance between molecules in microns
"""
function calc_r(mol1::DiffusingMolecule{<:Emitter3D}, mol2::DiffusingMolecule{<:Emitter3D})
    sqrt((mol1.x - mol2.x)^2 + (mol1.y - mol2.y)^2 + (mol1.z - mol2.z)^2)
end

"""
    calc_ϕ(mol1::DiffusingMolecule, mol2::DiffusingMolecule)

Calculate azimuthal angle between molecules.

# Arguments
- `mol1::DiffusingMolecule`: First molecule
- `mol2::DiffusingMolecule`: Second molecule

# Returns
- `Float64`: Azimuthal angle in radians
"""
function calc_ϕ(mol1::DiffusingMolecule, mol2::DiffusingMolecule)
    atan(mol2.y - mol1.y, mol2.x - mol1.x)
end

"""
    calc_θ(mol1::DiffusingMolecule{<:Emitter2D}, mol2::DiffusingMolecule{<:Emitter2D})

Calculate polar angle between 2D molecules (always π/2 for 2D).

# Arguments
- `mol1::DiffusingMolecule{<:Emitter2D}`: First molecule
- `mol2::DiffusingMolecule{<:Emitter2D}`: Second molecule

# Returns
- `Float64`: Polar angle in radians (always π/2 for 2D)
"""
function calc_θ(mol1::DiffusingMolecule{<:Emitter2D}, mol2::DiffusingMolecule{<:Emitter2D})
    π/2  # Always in xy plane for 2D
end

"""
    calc_θ(mol1::DiffusingMolecule{<:Emitter3D}, mol2::DiffusingMolecule{<:Emitter3D})

Calculate polar angle between 3D molecules.

# Arguments
- `mol1::DiffusingMolecule{<:Emitter3D}`: First molecule
- `mol2::DiffusingMolecule{<:Emitter3D}`: Second molecule

# Returns
- `Float64`: Polar angle in radians
"""
function calc_θ(mol1::DiffusingMolecule{<:Emitter3D}, mol2::DiffusingMolecule{<:Emitter3D})
    r = calc_r(mol1, mol2)
    acos((mol2.z - mol1.z) / r)
end

# Core molecular state changes
"""
    dimerize!(mol1::DiffusingMolecule, mol2::DiffusingMolecule, distance::Real)

Form a dimer between two molecules.

# Arguments
- `mol1::DiffusingMolecule`: First molecule
- `mol2::DiffusingMolecule`: Second molecule
- `distance::Real`: Distance between molecules in the dimer in microns

# Returns
- `Nothing`

# Example
```julia
dimerize!(molecule1, molecule2, 0.05)  # Form dimer with 50nm separation
```

# Note
Molecules are repositioned to maintain the specified distance.
The center of mass position is preserved.
"""
function dimerize!(mol1::DiffusingMolecule, mol2::DiffusingMolecule, distance::Real)
    # Input validation
    if mol1.state == 2 || mol2.state == 2
        return nothing  # Already part of dimers
    end
    
    if distance <= 0
        throw(ArgumentError("Distance must be positive"))
    end
    
    # Update states
    mol1.state = 2
    mol2.state = 2
    mol1.link = mol2.id
    mol2.link = mol1.id
    
    # Calculate center of mass
    com_x = (mol1.x + mol2.x) / 2
    com_y = (mol1.y + mol2.y) / 2
    
    # Calculate orientation
    ϕ = calc_ϕ(mol1, mol2)
    r = distance / 2
    
    # Set new positions based on COM and orientation
    dx = r * cos(ϕ)
    dy = r * sin(ϕ)
    
    # Update positions
    mol1.x = com_x - dx
    mol1.y = com_y - dy
    mol2.x = com_x + dx
    mol2.y = com_y + dy
    
    return nothing
end

"""
    monomerize!(mol::DiffusingMolecule, system::DiffusingMoleculeSystem)

Convert a dimer back to monomers.

# Arguments
- `mol::DiffusingMolecule`: Molecule that is part of a dimer
- `system::DiffusingMoleculeSystem`: System containing all molecules

# Returns
- `Nothing`

# Example
```julia
monomerize!(molecule, system)  # Break dimer into monomers
```

# Note
Both molecules in the dimer are changed to monomer state.
Positions are not changed.
"""
function monomerize!(mol::DiffusingMolecule, system::DiffusingMoleculeSystem)
    if mol.state != 2 || isnothing(mol.link)
        return nothing
    end
    
    # Find linked molecule
    linked_idx = findfirst(m -> m.id == mol.link, system.molecules)
    isnothing(linked_idx) && return nothing
    
    linked_mol = system.molecules[linked_idx]
    
    # Update states
    mol.state = 1
    linked_mol.state = 1
    mol.link = nothing
    linked_mol.link = nothing
    
    return nothing
end