
abstract type AbstractOligomer end


"""
    struct Monomer <: AbstractOligomer

A struct representing a single monomer in an oligomer.

# Fields
- `x::Float64`: the x-coordinate of the monomer's position
- `y::Float64`: the y-coordinate of the monomer's position
- `z::Float64`: the z-coordinate of the monomer's position
- `state::Int64`: the state of the monomer (1 for monomer, 2 for dimer)
- `link::Union{Monomer,Nothing}`: the monomer that this monomer is linked to, if any
- `updated::Bool`: a flag indicating whether this monomer has been updated in the current iteration
"""
mutable struct Monomer <: AbstractOligomer
    x::Float64
    y::Float64
    z::Float64
    state::Int64
    link::Union{Monomer,Nothing}
    updated::Bool
end

"""
    MoleculeFrame(frame::Int64, molecules::Vector{<:AbstractOligomer})

A struct representing a time frame of a simulation of a system of molecules.

# Fields
- `frame::Int64`: The frame number of the simulation.
- `molecules::Vector{<:AbstractOligomer}`: A vector of `AbstractOligomer` objects representing the molecules in the system at the given frame.
"""
struct MoleculeFrame
    framenum::Int64
    molecules::Vector{<:AbstractOligomer}
end

function MoleculeFrame(framenum::Int64, nmolecules::Int64)
    molecules = [Monomer(0.0, 0.0, 0.0, 1, nothing, false) for i in 1:nmolecules]
    return MoleculeFrame(framenum, molecules)
end

"""
    MoleculeFrame(frame::Int64, molecules::Vector{<:AbstractOligomer})

A struct representing a frame of a simulation of a system of molecules.

# Fields
- `frame::Int64`: The frame number of the simulation.
- `molecules::Vector{<:AbstractOligomer}`: A vector of `AbstractOligomer` objects representing the molecules in the system at the given frame.
"""
struct MoleculeHistory
    dt::Float64
    frames::Vector{MoleculeFrame}
end

function MoleculeHistory(dt::Float64, nframes::Int64, nmolecules::Int64)
    frames = [MoleculeFrame(i, [Monomer(0.0, 0.0, 0.0, 1, nothing, false) for j in 1:nmolecules]) for i in 1:nframes]
    return MoleculeHistory(dt, frames)
end

