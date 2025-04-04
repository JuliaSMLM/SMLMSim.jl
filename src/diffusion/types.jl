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