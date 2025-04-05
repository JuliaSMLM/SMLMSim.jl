# Define specialized emitter types for diffusion simulation that
# directly inherit from SMLMData.AbstractEmitter

"""
    DiffusingEmitter2D{T<:AbstractFloat} <: AbstractEmitter

A 2D emitter type for diffusion simulations that contains both spatial
and temporal information, plus molecular state information.

# Fields
- `x::T`: x-coordinate in microns
- `y::T`: y-coordinate in microns
- `photons::T`: number of photons emitted
- `timestamp::T`: actual simulation time in seconds
- `frame::Int`: camera frame number based on framerate and exposure
- `dataset::Int`: dataset identifier
- `track_id::Int`: unique molecule identifier for tracking
- `state::Int`: molecular state (1=monomer, 2=dimer)
- `link_id::Union{Int,Nothing}`: ID of linked molecule (for dimers), or nothing for monomers
"""
struct DiffusingEmitter2D{T<:AbstractFloat} <: AbstractEmitter
    # Core spatial and intensity properties required by integrate_pixels
    x::T
    y::T
    photons::T
    
    # Temporal properties
    timestamp::T  # Actual simulation time in seconds
    frame::Int    # Frame number based on camera integration
    
    # Bookkeeping properties
    dataset::Int
    track_id::Int
    
    # Diffusion-specific properties
    state::Int              # 1=monomer, 2=dimer
    link_id::Union{Int,Nothing}   # ID of linked molecule (for dimers)
end

"""
    DiffusingEmitter3D{T<:AbstractFloat} <: AbstractEmitter

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
- `track_id::Int`: unique molecule identifier for tracking
- `state::Int`: molecular state (1=monomer, 2=dimer)
- `link_id::Union{Int,Nothing}`: ID of linked molecule (for dimers), or nothing for monomers
"""
struct DiffusingEmitter3D{T<:AbstractFloat} <: AbstractEmitter
    # Core spatial and intensity properties
    x::T
    y::T
    z::T
    photons::T
    
    # Same temporal and diffusion properties
    timestamp::T
    frame::Int
    dataset::Int
    track_id::Int
    state::Int
    link_id::Union{Int,Nothing}
end

# Display methods
function Base.show(io::IO, e::DiffusingEmitter2D{T}) where T
    state_str = e.state == 1 ? "monomer" : "dimer"
    print(io, "DiffusingEmitter2D{$T}($(e.x), $(e.y) μm, t=$(e.timestamp)s, frame=$(e.frame), $state_str)")
end

function Base.show(io::IO, ::MIME"text/plain", e::DiffusingEmitter2D{T}) where T
    state_str = e.state == 1 ? "monomer" : "dimer"
    link_str = isnothing(e.link_id) ? "unlinked" : "linked to $(e.link_id)"
    
    println(io, "DiffusingEmitter2D{$T}:")
    println(io, "  Position: ($(e.x), $(e.y)) μm")
    println(io, "  Photons: $(e.photons)")
    println(io, "  Time: $(e.timestamp) s, frame: $(e.frame)")
    println(io, "  State: $state_str")
    println(io, "  Track ID: $(e.track_id)")
    print(io, "  Link: $link_str")
end

function Base.show(io::IO, e::DiffusingEmitter3D{T}) where T
    state_str = e.state == 1 ? "monomer" : "dimer"
    print(io, "DiffusingEmitter3D{$T}($(e.x), $(e.y), $(e.z) μm, t=$(e.timestamp)s, frame=$(e.frame), $state_str)")
end

function Base.show(io::IO, ::MIME"text/plain", e::DiffusingEmitter3D{T}) where T
    state_str = e.state == 1 ? "monomer" : "dimer"
    link_str = isnothing(e.link_id) ? "unlinked" : "linked to $(e.link_id)"
    
    println(io, "DiffusingEmitter3D{$T}:")
    println(io, "  Position: ($(e.x), $(e.y), $(e.z)) μm")
    println(io, "  Photons: $(e.photons)")
    println(io, "  Time: $(e.timestamp) s, frame: $(e.frame)")
    println(io, "  State: $state_str")
    println(io, "  Track ID: $(e.track_id)")
    print(io, "  Link: $link_str")
end
