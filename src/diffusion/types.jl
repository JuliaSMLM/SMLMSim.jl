# Define specialized emitter types for diffusion simulation that
# directly inherit from SMLMData.AbstractEmitter

"""
    AbstractDiffusingEmitter <: AbstractEmitter

Abstract type for all diffusing emitters to enable dispatch-based operations.
This provides a common parent for 2D and 3D diffusing emitters.
"""
abstract type AbstractDiffusingEmitter <: AbstractEmitter end

"""
    DiffusingEmitter2D{T<:AbstractFloat} <: AbstractDiffusingEmitter

A 2D emitter type for diffusion simulations that contains both spatial
and temporal information, plus molecular state information.

# Fields
- `x::T`: x-coordinate in microns
- `y::T`: y-coordinate in microns
- `photons::T`: number of photons emitted
- `timestamp::T`: actual simulation time in seconds
- `frame::Int`: camera frame number based on framerate and exposure
- `dataset::Int`: dataset identifier
- `track_id::Int`: unique molecule identifier 
- `state::Symbol`: molecular state (:monomer or :dimer)
- `partner_id::Union{Int,Nothing}`: ID of linked molecule (for dimers), or nothing for monomers
"""
struct DiffusingEmitter2D{T<:AbstractFloat} <: AbstractDiffusingEmitter
    # Core spatial and intensity properties required by integrate_pixels
    x::T
    y::T
    photons::T
    
    # Temporal properties
    timestamp::T  # Actual simulation time in seconds
    frame::Int    # Frame number based on camera integration
    
    # Bookkeeping properties
    dataset::Int
    track_id::Int  # Changed from id to track_id to be consistent with other emitters
    
    # Diffusion-specific properties
    state::Symbol              # :monomer or :dimer
    partner_id::Union{Int,Nothing}   # ID of linked molecule (for dimers)
end

"""
    DiffusingEmitter3D{T<:AbstractFloat} <: AbstractDiffusingEmitter

A 3D emitter type for diffusion simulations that contains both spatial
and temporal information, plus molecular state information.

# Fields
- `x::T`: x-coordinate in microns
- `y::T`: y-coordinate in microns
- `z::T`: z-coordinate in microns
- `photons::T`: number of photons emitted
- `timestamp::T`: actual simulation time in seconds
- `frame::Int`: camera frame number based on framerate and exposure
- `dataset::Int`: dataset identifier
- `track_id::Int`: unique molecule identifier
- `state::Symbol`: molecular state (:monomer or :dimer)
- `partner_id::Union{Int,Nothing}`: ID of linked molecule (for dimers), or nothing for monomers
"""
struct DiffusingEmitter3D{T<:AbstractFloat} <: AbstractDiffusingEmitter
    # Core spatial and intensity properties
    x::T
    y::T
    z::T
    photons::T
    
    # Same temporal and diffusion properties
    timestamp::T
    frame::Int
    dataset::Int
    track_id::Int  # Changed from id to track_id to be consistent with other emitters
    state::Symbol
    partner_id::Union{Int,Nothing}
end

# Display methods
function Base.show(io::IO, e::DiffusingEmitter2D{T}) where T
    print(io, "DiffusingEmitter2D{$T}($(e.x), $(e.y) μm, t=$(e.timestamp)s, frame=$(e.frame), $(e.state))")
end

function Base.show(io::IO, ::MIME"text/plain", e::DiffusingEmitter2D{T}) where T
    link_str = isnothing(e.partner_id) ? "unlinked" : "linked to $(e.partner_id)"
    
    println(io, "DiffusingEmitter2D{$T}:")
    println(io, "  Position: ($(e.x), $(e.y)) μm")
    println(io, "  Photons: $(e.photons)")
    println(io, "  Time: $(e.timestamp) s, frame: $(e.frame)")
    println(io, "  State: $(e.state)")
    println(io, "  ID: $(e.track_id)")
    print(io, "  Link: $link_str")
end

function Base.show(io::IO, e::DiffusingEmitter3D{T}) where T
    print(io, "DiffusingEmitter3D{$T}($(e.x), $(e.y), $(e.z) μm, t=$(e.timestamp)s, frame=$(e.frame), $(e.state))")
end

function Base.show(io::IO, ::MIME"text/plain", e::DiffusingEmitter3D{T}) where T
    link_str = isnothing(e.partner_id) ? "unlinked" : "linked to $(e.partner_id)"
    
    println(io, "DiffusingEmitter3D{$T}:")
    println(io, "  Position: ($(e.x), $(e.y), $(e.z)) μm")
    println(io, "  Photons: $(e.photons)")
    println(io, "  Time: $(e.timestamp) s, frame: $(e.frame)")
    println(io, "  State: $(e.state)")
    println(io, "  ID: $(e.track_id)")
    print(io, "  Link: $link_str")
end

# Add copy and deepcopy methods for diffusing emitter types
"""
    Base.copy(e::DiffusingEmitter2D)

Create a copy of a 2D diffusing emitter.
"""
function Base.copy(e::DiffusingEmitter2D{T}) where T <: AbstractFloat
    DiffusingEmitter2D{T}(
        e.x, e.y,           # Position
        e.photons,          # Photons
        e.timestamp,        # Timestamp
        e.frame,            # Frame
        e.dataset,          # Dataset
        e.track_id,         # ID (now track_id)
        e.state,            # State
        e.partner_id        # Partner ID
    )
end

"""
    Base.copy(e::DiffusingEmitter3D)

Create a copy of a 3D diffusing emitter.
"""
function Base.copy(e::DiffusingEmitter3D{T}) where T <: AbstractFloat
    DiffusingEmitter3D{T}(
        e.x, e.y, e.z,      # Position
        e.photons,          # Photons
        e.timestamp,        # Timestamp
        e.frame,            # Frame
        e.dataset,          # Dataset
        e.track_id,         # ID (now track_id)
        e.state,            # State
        e.partner_id        # Partner ID
    )
end

# Add deepcopy methods (for immutable structs, copy and deepcopy can be the same)
Base.deepcopy(e::DiffusingEmitter2D) = copy(e)
Base.deepcopy(e::DiffusingEmitter3D) = copy(e)
