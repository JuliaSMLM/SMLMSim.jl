# src/types.jl

"""
    DiffusingMolecule{E<:AbstractEmitter}

A molecule that can diffuse and form dimers, containing an SMLM emitter.

# Fields
- `emitter::E`: underlying SMLM emitter (2D/3D, fitted/unfitted)
- `state::Int`: oligomerization state (1=monomer, 2=dimer)
- `id::Int`: unique identifier for tracking
- `link::Union{Nothing,Int}`: id of linked molecule (for dimers)
- `updated::Bool`: flag for position updates
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
"""Calculate Euclidean distance between two molecules"""
function calc_r(mol1::DiffusingMolecule{<:Emitter2D}, mol2::DiffusingMolecule{<:Emitter2D})
    sqrt((mol1.x - mol2.x)^2 + (mol1.y - mol2.y)^2)
end

function calc_r(mol1::DiffusingMolecule{<:Emitter3D}, mol2::DiffusingMolecule{<:Emitter3D})
    sqrt((mol1.x - mol2.x)^2 + (mol1.y - mol2.y)^2 + (mol1.z - mol2.z)^2)
end

"""Calculate azimuthal angle between molecules"""
function calc_ϕ(mol1::DiffusingMolecule, mol2::DiffusingMolecule)
    atan(mol2.y - mol1.y, mol2.x - mol1.x)
end

"""Calculate polar angle between molecules"""
function calc_θ(mol1::DiffusingMolecule{<:Emitter2D}, mol2::DiffusingMolecule{<:Emitter2D})
    π/2  # Always in xy plane for 2D
end

function calc_θ(mol1::DiffusingMolecule{<:Emitter3D}, mol2::DiffusingMolecule{<:Emitter3D})
    r = calc_r(mol1, mol2)
    acos((mol2.z - mol1.z) / r)
end

# Core molecular state changes
"""
    dimerize!(mol1::DiffusingMolecule, mol2::DiffusingMolecule, distance::Real)

Form a dimer between two molecules.
"""
function dimerize!(mol1::DiffusingMolecule, mol2::DiffusingMolecule, distance::Real)
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