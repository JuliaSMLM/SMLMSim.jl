"""
    Helpers for diffusion simulation.

This file contains helper functions for the diffusion simulation module,
including distance calculations and state management using dispatch-based operations.
"""

# Helper functions for geometric calculations
"""
    distance(e1, e2)

Calculate Euclidean distance between two emitters.
Generic implementation with multiple dispatch.

# Arguments
- `e1`: First emitter
- `e2`: Second emitter

# Returns
- `Float64`: Distance between emitters in microns
"""
# Generic fallback method
function distance(e1, e2)
    error("No distance method implemented for $(typeof(e1)) and $(typeof(e2))")
end

# Specialized for 2D emitters
function distance(e1::DiffusingEmitter2D{T}, e2::DiffusingEmitter2D{T}) where T <: AbstractFloat
    sqrt((e1.x - e2.x)^2 + (e1.y - e2.y)^2)
end

# Specialized for 3D emitters
function distance(e1::DiffusingEmitter3D{T}, e2::DiffusingEmitter3D{T}) where T <: AbstractFloat
    sqrt((e1.x - e2.x)^2 + (e1.y - e2.y)^2 + (e1.z - e2.z)^2)
end

# Mixed dimensions (fallback)
function distance(e1::AbstractDiffusingEmitter, e2::AbstractDiffusingEmitter)
    error("Cannot calculate distance between different emitter dimensions")
end

"""
    angle(e1, e2)

Calculate angle between emitters.
Generic implementation with multiple dispatch.

# Arguments
- `e1`: First emitter
- `e2`: Second emitter

# Returns
- Angle representation appropriate for the emitter dimensions
"""
# Generic fallback method
function angle(e1, e2)
    error("No angle method implemented for $(typeof(e1)) and $(typeof(e2))")
end

# Specialized for 2D emitters - returns azimuthal angle
function angle(e1::DiffusingEmitter2D{T}, e2::DiffusingEmitter2D{T}) where T <: AbstractFloat
    atan(e2.y - e1.y, e2.x - e1.x)
end

# Specialized for 3D emitters - returns (azimuthal, polar) angles
function angle(e1::DiffusingEmitter3D{T}, e2::DiffusingEmitter3D{T}) where T <: AbstractFloat
    # Azimuthal angle (ϕ)
    ϕ = atan(e2.y - e1.y, e2.x - e1.x)
    
    # Polar angle (θ)
    r = distance(e1, e2)
    θ = acos((e2.z - e1.z) / r)
    
    return (ϕ, θ)
end

# More generic approach with coordinate tuples
function distance(p1::NTuple{N,T}, p2::NTuple{N,T}) where {N,T<:AbstractFloat}
    sqrt(sum((p1[i] - p2[i])^2 for i in 1:N))
end

# Extract coordinates as tuples
function coordinates(e::DiffusingEmitter2D{T}) where T <: AbstractFloat
    (e.x, e.y)
end

function coordinates(e::DiffusingEmitter3D{T}) where T <: AbstractFloat
    (e.x, e.y, e.z)
end

# State management functions
"""
    can_dimerize(e1::AbstractDiffusingEmitter, e2::AbstractDiffusingEmitter, r_react::Float64)

Check if two emitters can form a dimer.

# Arguments
- `e1::AbstractDiffusingEmitter`: First emitter
- `e2::AbstractDiffusingEmitter`: Second emitter
- `r_react::Float64`: Reaction radius in microns

# Returns
- `Bool`: True if emitters can form a dimer
"""
function can_dimerize(e1::AbstractDiffusingEmitter, e2::AbstractDiffusingEmitter, r_react::Float64)
    e1.state == :monomer && 
    e2.state == :monomer && 
    distance(e1, e2) < r_react
end

"""
    dimerize(e1::DiffusingEmitter2D, e2::DiffusingEmitter2D, d_dimer::Float64)

Create two new emitters in dimer state from two monomers.

# Arguments
- `e1::DiffusingEmitter2D`: First emitter
- `e2::DiffusingEmitter2D`: Second emitter
- `d_dimer::Float64`: Dimer separation distance in microns

# Returns
- `Tuple{DiffusingEmitter2D, DiffusingEmitter2D}`: Two new emitters in dimer state
"""
function dimerize(e1::DiffusingEmitter2D{T}, e2::DiffusingEmitter2D{T}, d_dimer::Float64) where T <: AbstractFloat
    # Calculate center of mass
    com_x = (e1.x + e2.x) / 2
    com_y = (e1.y + e2.y) / 2
    
    # Calculate orientation
    ϕ = angle(e1, e2)
    r = d_dimer / 2
    
    # Calculate new positions
    dx = r * cos(ϕ)
    dy = r * sin(ϕ)
    
    # Create new dimer emitters
    d1 = DiffusingEmitter2D{T}(
        com_x - dx, com_y - dy,  # Position
        e1.photons,              # Photons
        e1.timestamp,            # Timestamp
        e1.frame,                # Frame
        e1.dataset,              # Dataset
        e1.id,                   # ID
        :dimer,                  # State
        e2.id                    # Partner ID
    )
    
    d2 = DiffusingEmitter2D{T}(
        com_x + dx, com_y + dy,  # Position
        e2.photons,              # Photons
        e2.timestamp,            # Timestamp
        e2.frame,                # Frame
        e2.dataset,              # Dataset
        e2.id,                   # ID
        :dimer,                  # State
        e1.id                    # Partner ID
    )
    
    return (d1, d2)
end

"""
    dimerize(e1::DiffusingEmitter3D, e2::DiffusingEmitter3D, d_dimer::Float64)

Create two new emitters in dimer state from two monomers in 3D.

# Arguments
- `e1::DiffusingEmitter3D`: First emitter
- `e2::DiffusingEmitter3D`: Second emitter
- `d_dimer::Float64`: Dimer separation distance in microns

# Returns
- `Tuple{DiffusingEmitter3D, DiffusingEmitter3D}`: Two new emitters in dimer state
"""
function dimerize(e1::DiffusingEmitter3D{T}, e2::DiffusingEmitter3D{T}, d_dimer::Float64) where T <: AbstractFloat
    # Calculate center of mass
    com_x = (e1.x + e2.x) / 2
    com_y = (e1.y + e2.y) / 2
    com_z = (e1.z + e2.z) / 2
    
    # Calculate orientation
    ϕ, θ = angle(e1, e2)
    r = d_dimer / 2
    
    # Calculate new positions
    dx = r * sin(θ) * cos(ϕ)
    dy = r * sin(θ) * sin(ϕ)
    dz = r * cos(θ)
    
    # Create new dimer emitters
    d1 = DiffusingEmitter3D{T}(
        com_x - dx, com_y - dy, com_z - dz,  # Position
        e1.photons,                          # Photons
        e1.timestamp,                        # Timestamp
        e1.frame,                            # Frame
        e1.dataset,                          # Dataset
        e1.id,                               # ID
        :dimer,                              # State
        e2.id                                # Partner ID
    )
    
    d2 = DiffusingEmitter3D{T}(
        com_x + dx, com_y + dy, com_z + dz,  # Position
        e2.photons,                          # Photons
        e2.timestamp,                        # Timestamp
        e2.frame,                            # Frame
        e2.dataset,                          # Dataset
        e2.id,                               # ID
        :dimer,                              # State
        e1.id                                # Partner ID
    )
    
    return (d1, d2)
end

"""
    should_dissociate(e::AbstractDiffusingEmitter, k_off::Float64, dt::Float64)

Check if a dimer should dissociate based on stochastic rate.

# Arguments
- `e::AbstractDiffusingEmitter`: Emitter to check
- `k_off::Float64`: Dissociation rate (s⁻¹)
- `dt::Float64`: Time step (s)

# Returns
- `Bool`: True if dimer should dissociate
"""
function should_dissociate(e::AbstractDiffusingEmitter, k_off::Float64, dt::Float64)
    e.state == :dimer && rand() < k_off * dt
end

"""
    dissociate(e::DiffusingEmitter2D, emitters::Vector{<:AbstractDiffusingEmitter})

Create two new monomers from a dimer.

# Arguments
- `e::DiffusingEmitter2D`: Emitter part of a dimer
- `emitters::Vector{<:AbstractDiffusingEmitter}`: All emitters in the system

# Returns
- `Tuple{DiffusingEmitter2D, DiffusingEmitter2D}`: Two new emitters in monomer state
"""
function dissociate(e::DiffusingEmitter2D{T}, emitters::Vector{<:AbstractDiffusingEmitter}) where T <: AbstractFloat
    # Find partner
    if isnothing(e.partner_id)
        error("Emitter is not part of a dimer")
    end
    
    partner_idx = findfirst(em -> em.id == e.partner_id, emitters)
    if isnothing(partner_idx)
        error("Partner emitter not found")
    end
    
    partner = emitters[partner_idx]
    
    # Create new monomer emitters at same positions
    m1 = DiffusingEmitter2D{T}(
        e.x, e.y,           # Position
        e.photons,          # Photons
        e.timestamp,        # Timestamp
        e.frame,            # Frame
        e.dataset,          # Dataset
        e.id,               # ID
        :monomer,           # State
        nothing             # Partner ID
    )
    
    m2 = DiffusingEmitter2D{T}(
        partner.x, partner.y,  # Position
        partner.photons,       # Photons
        partner.timestamp,     # Timestamp
        partner.frame,         # Frame
        partner.dataset,       # Dataset
        partner.id,            # ID
        :monomer,              # State
        nothing                # Partner ID
    )
    
    return (m1, m2)
end

"""
    dissociate(e::DiffusingEmitter3D, emitters::Vector{<:AbstractDiffusingEmitter})

Create two new monomers from a 3D dimer.

# Arguments
- `e::DiffusingEmitter3D`: Emitter part of a dimer
- `emitters::Vector{<:AbstractDiffusingEmitter}`: All emitters in the system

# Returns
- `Tuple{DiffusingEmitter3D, DiffusingEmitter3D}`: Two new emitters in monomer state
"""
function dissociate(e::DiffusingEmitter3D{T}, emitters::Vector{<:AbstractDiffusingEmitter}) where T <: AbstractFloat
    # Find partner
    if isnothing(e.partner_id)
        error("Emitter is not part of a dimer")
    end
    
    partner_idx = findfirst(em -> em.id == e.partner_id, emitters)
    if isnothing(partner_idx)
        error("Partner emitter not found")
    end
    
    partner = emitters[partner_idx]
    
    # Create new monomer emitters at same positions
    m1 = DiffusingEmitter3D{T}(
        e.x, e.y, e.z,      # Position
        e.photons,          # Photons
        e.timestamp,        # Timestamp
        e.frame,            # Frame
        e.dataset,          # Dataset
        e.id,               # ID
        :monomer,           # State
        nothing             # Partner ID
    )
    
    m2 = DiffusingEmitter3D{T}(
        partner.x, partner.y, partner.z,  # Position
        partner.photons,                  # Photons
        partner.timestamp,                # Timestamp
        partner.frame,                    # Frame
        partner.dataset,                  # Dataset
        partner.id,                       # ID
        :monomer,                         # State
        nothing                           # Partner ID
    )
    
    return (m1, m2)
end

# Diffusion functions
"""
    diffuse(e::DiffusingEmitter2D, diff_coef::Float64, dt::Float64)

Create a new emitter with updated position based on Brownian motion.

# Arguments
- `e::DiffusingEmitter2D`: Emitter to update
- `diff_coef::Float64`: Diffusion coefficient (μm²/s)
- `dt::Float64`: Time step (s)

# Returns
- `DiffusingEmitter2D`: New emitter with updated position
"""
function diffuse(e::DiffusingEmitter2D{T}, diff_coef::Float64, dt::Float64) where T <: AbstractFloat
    σ = sqrt(2 * diff_coef * dt)
    
    # Apply Brownian motion
    new_x = e.x + rand(Normal(0, σ))
    new_y = e.y + rand(Normal(0, σ))
    
    # Create new emitter with updated position
    DiffusingEmitter2D{T}(
        new_x, new_y,       # Updated position
        e.photons,          # Photons
        e.timestamp + dt,   # Updated timestamp
        e.frame,            # Frame
        e.dataset,          # Dataset
        e.id,               # ID
        e.state,            # State
        e.partner_id        # Partner ID
    )
end

"""
    diffuse(e::DiffusingEmitter3D, diff_coef::Float64, dt::Float64)

Create a new 3D emitter with updated position based on Brownian motion.

# Arguments
- `e::DiffusingEmitter3D`: Emitter to update
- `diff_coef::Float64`: Diffusion coefficient (μm²/s)
- `dt::Float64`: Time step (s)

# Returns
- `DiffusingEmitter3D`: New emitter with updated position
"""
function diffuse(e::DiffusingEmitter3D{T}, diff_coef::Float64, dt::Float64) where T <: AbstractFloat
    σ = sqrt(2 * diff_coef * dt)
    
    # Apply Brownian motion
    new_x = e.x + rand(Normal(0, σ))
    new_y = e.y + rand(Normal(0, σ))
    new_z = e.z + rand(Normal(0, σ))
    
    # Create new emitter with updated position
    DiffusingEmitter3D{T}(
        new_x, new_y, new_z,  # Updated position
        e.photons,            # Photons
        e.timestamp + dt,     # Updated timestamp
        e.frame,              # Frame
        e.dataset,            # Dataset
        e.id,                 # ID
        e.state,              # State
        e.partner_id          # Partner ID
    )
end

"""
    diffuse_dimer(e1::DiffusingEmitter2D, e2::DiffusingEmitter2D, diff_trans::Float64, diff_rot::Float64, d_dimer::Float64, dt::Float64)

Diffuse a dimer with both translational and rotational components.

# Arguments
- `e1::DiffusingEmitter2D`: First emitter in dimer
- `e2::DiffusingEmitter2D`: Second emitter in dimer
- `diff_trans::Float64`: Translational diffusion coefficient (μm²/s)
- `diff_rot::Float64`: Rotational diffusion coefficient (rad²/s)
- `d_dimer::Float64`: Dimer separation distance (μm)
- `dt::Float64`: Time step (s)

# Returns
- `Tuple{DiffusingEmitter2D, DiffusingEmitter2D}`: Two new emitters with updated positions
"""
function diffuse_dimer(e1::DiffusingEmitter2D{T}, e2::DiffusingEmitter2D{T}, diff_trans::Float64, diff_rot::Float64, d_dimer::Float64, dt::Float64) where T <: AbstractFloat
    # Translational diffusion
    σ_trans = sqrt(2 * diff_trans * dt)
    dx = rand(Normal(0, σ_trans))
    dy = rand(Normal(0, σ_trans))
    
    # Calculate center of mass
    com_x = (e1.x + e2.x) / 2
    com_y = (e1.y + e2.y) / 2
    
    # Move center of mass
    com_x += dx
    com_y += dy
    
    # Rotational diffusion
    σ_rot = sqrt(2 * diff_rot * dt)
    ϕ = angle(e1, e2)
    ϕ += rand(Normal(0, σ_rot))
    
    # Calculate new positions
    r = d_dimer / 2
    dx = r * cos(ϕ)
    dy = r * sin(ϕ)
    
    # Create new emitters with updated positions
    d1 = DiffusingEmitter2D{T}(
        com_x - dx, com_y - dy,  # Updated position
        e1.photons,              # Photons
        e1.timestamp + dt,       # Updated timestamp
        e1.frame,                # Frame
        e1.dataset,              # Dataset
        e1.id,                   # ID
        e1.state,                # State
        e1.partner_id            # Partner ID
    )
    
    d2 = DiffusingEmitter2D{T}(
        com_x + dx, com_y + dy,  # Updated position
        e2.photons,              # Photons
        e2.timestamp + dt,       # Updated timestamp
        e2.frame,                # Frame
        e2.dataset,              # Dataset
        e2.id,                   # ID
        e2.state,                # State
        e2.partner_id            # Partner ID
    )
    
    return (d1, d2)
end

"""
    diffuse_dimer(e1::DiffusingEmitter3D, e2::DiffusingEmitter3D, diff_trans::Float64, diff_rot::Float64, d_dimer::Float64, dt::Float64)

Diffuse a 3D dimer with both translational and rotational components.

# Arguments
- `e1::DiffusingEmitter3D`: First emitter in dimer
- `e2::DiffusingEmitter3D`: Second emitter in dimer
- `diff_trans::Float64`: Translational diffusion coefficient (μm²/s)
- `diff_rot::Float64`: Rotational diffusion coefficient (rad²/s)
- `d_dimer::Float64`: Dimer separation distance (μm)
- `dt::Float64`: Time step (s)

# Returns
- `Tuple{DiffusingEmitter3D, DiffusingEmitter3D}`: Two new emitters with updated positions
"""
function diffuse_dimer(e1::DiffusingEmitter3D{T}, e2::DiffusingEmitter3D{T}, diff_trans::Float64, diff_rot::Float64, d_dimer::Float64, dt::Float64) where T <: AbstractFloat
    # Translational diffusion
    σ_trans = sqrt(2 * diff_trans * dt)
    dx = rand(Normal(0, σ_trans))
    dy = rand(Normal(0, σ_trans))
    dz = rand(Normal(0, σ_trans))
    
    # Calculate center of mass
    com_x = (e1.x + e2.x) / 2
    com_y = (e1.y + e2.y) / 2
    com_z = (e1.z + e2.z) / 2
    
    # Move center of mass
    com_x += dx
    com_y += dy
    com_z += dz
    
    # Rotational diffusion
    σ_rot = sqrt(2 * diff_rot * dt)
    ϕ, θ = angle(e1, e2)
    
    # Apply random rotation (simplified model)
    ϕ += rand(Normal(0, σ_rot))
    θ += rand(Normal(0, σ_rot))
    
    # Calculate new positions
    r = d_dimer / 2
    dx = r * sin(θ) * cos(ϕ)
    dy = r * sin(θ) * sin(ϕ)
    dz = r * cos(θ)
    
    # Create new emitters with updated positions
    d1 = DiffusingEmitter3D{T}(
        com_x - dx, com_y - dy, com_z - dz,  # Updated position
        e1.photons,                          # Photons
        e1.timestamp + dt,                   # Updated timestamp
        e1.frame,                            # Frame
        e1.dataset,                          # Dataset
        e1.id,                               # ID
        e1.state,                            # State
        e1.partner_id                        # Partner ID
    )
    
    d2 = DiffusingEmitter3D{T}(
        com_x + dx, com_y + dy, com_z + dz,  # Updated position
        e2.photons,                          # Photons
        e2.timestamp + dt,                   # Updated timestamp
        e2.frame,                            # Frame
        e2.dataset,                          # Dataset
        e2.id,                               # ID
        e2.state,                            # State
        e2.partner_id                        # Partner ID
    )
    
    return (d1, d2)
end

"""
    apply_boundary(e::AbstractDiffusingEmitter, box_size::Float64, boundary::String)

Apply boundary conditions to an emitter (generic fallback).

# Arguments
- `e::AbstractDiffusingEmitter`: Emitter to apply boundary to
- `box_size::Float64`: Simulation box size in microns
- `boundary::String`: Boundary condition type ("periodic" or "reflecting")

# Returns
- `AbstractDiffusingEmitter`: Emitter with boundary conditions applied
"""
function apply_boundary(e::AbstractDiffusingEmitter, box_size::Float64, boundary::String)
    error("No boundary method implemented for $(typeof(e))")
end

"""
    apply_boundary(e::DiffusingEmitter2D, box_size::Float64, boundary::String)

Apply boundary conditions to a 2D emitter.

# Arguments
- `e::DiffusingEmitter2D`: Emitter to apply boundary to
- `box_size::Float64`: Simulation box size in microns
- `boundary::String`: Boundary condition type ("periodic" or "reflecting")

# Returns
- `DiffusingEmitter2D`: New emitter with position constrained to the box
"""
function apply_boundary(e::DiffusingEmitter2D{T}, box_size::Float64, boundary::String) where T <: AbstractFloat
    new_x = e.x
    new_y = e.y
    
    if boundary == "periodic"
        new_x = mod(new_x, box_size)
        new_y = mod(new_y, box_size)
    else  # reflecting
        if new_x < 0
            new_x = -new_x
        elseif new_x > box_size
            new_x = 2box_size - new_x
        end
        
        if new_y < 0
            new_y = -new_y
        elseif new_y > box_size
            new_y = 2box_size - new_y
        end
    end
    
    # Only create a new emitter if the position changed
    if new_x != e.x || new_y != e.y
        return DiffusingEmitter2D{T}(
            new_x, new_y,      # Updated position
            e.photons,         # Photons
            e.timestamp,       # Timestamp
            e.frame,           # Frame
            e.dataset,         # Dataset
            e.id,              # ID
            e.state,           # State
            e.partner_id       # Partner ID
        )
    else
        return e
    end
end

"""
    apply_boundary(e::DiffusingEmitter3D, box_size::Float64, boundary::String)

Apply boundary conditions to a 3D emitter.

# Arguments
- `e::DiffusingEmitter3D`: Emitter to apply boundary to
- `box_size::Float64`: Simulation box size in microns
- `boundary::String`: Boundary condition type ("periodic" or "reflecting")

# Returns
- `DiffusingEmitter3D`: New emitter with position constrained to the box
"""
function apply_boundary(e::DiffusingEmitter3D{T}, box_size::Float64, boundary::String) where T <: AbstractFloat
    new_x = e.x
    new_y = e.y
    new_z = e.z
    
    if boundary == "periodic"
        new_x = mod(new_x, box_size)
        new_y = mod(new_y, box_size)
        new_z = mod(new_z, box_size)
    else  # reflecting
        if new_x < 0
            new_x = -new_x
        elseif new_x > box_size
            new_x = 2box_size - new_x
        end
        
        if new_y < 0
            new_y = -new_y
        elseif new_y > box_size
            new_y = 2box_size - new_y
        end
        
        if new_z < 0
            new_z = -new_z
        elseif new_z > box_size
            new_z = 2box_size - new_z
        end
    end
    
    # Only create a new emitter if the position changed
    if new_x != e.x || new_y != e.y || new_z != e.z
        return DiffusingEmitter3D{T}(
            new_x, new_y, new_z,  # Updated position
            e.photons,            # Photons
            e.timestamp,          # Timestamp
            e.frame,              # Frame
            e.dataset,            # Dataset
            e.id,                 # ID
            e.state,              # State
            e.partner_id          # Partner ID
        )
    else
        return e
    end
end

# SMLD conversion utilities
"""
    create_smld(emitters::Vector{<:AbstractDiffusingEmitter}, camera::AbstractCamera, params::DiffusionSMLMParams)

Convert a collection of diffusing emitters to a BasicSMLD object.

# Arguments
- `emitters::Vector{<:AbstractDiffusingEmitter}`: Collection of emitters from simulation
- `camera::AbstractCamera`: Camera model for imaging
- `params::DiffusionSMLMParams`: Simulation parameters

# Returns
- `BasicSMLD`: SMLD containing all emitters for further analysis or visualization
"""
function create_smld(emitters::Vector{<:AbstractDiffusingEmitter}, camera::AbstractCamera, params::DiffusionSMLMParams)
    # Determine max frame number
    max_frame = isempty(emitters) ? 0 : maximum(e -> e.frame, emitters)
    
    # Create metadata
    metadata = Dict{String,Any}(
        "simulation_type" => "diffusion",
        "simulation_parameters" => params,
        "camera_framerate" => params.camera_framerate,
        "camera_exposure" => params.camera_exposure
    )
    
    # Create SMLD object
    return BasicSMLD(
        emitters,
        camera,
        max_frame,       # nframes based on camera framerate
        1,               # ndatasets
        metadata
    )
end

"""
    get_frame(smld::BasicSMLD, frame_num::Int)

Extract emitters from a specific frame.

# Arguments
- `smld::BasicSMLD`: SMLD containing all emitters
- `frame_num::Int`: Frame number to extract

# Returns
- `BasicSMLD`: New SMLD containing only emitters from the specified frame
"""
function get_frame(smld::BasicSMLD, frame_num::Int)
    # Filter emitters by frame
    frame_emitters = filter(e -> e.frame == frame_num, smld.emitters)
    
    # Create new SMLD with same metadata
    return BasicSMLD(
        frame_emitters,
        smld.camera,
        1,  # Single frame
        smld.n_datasets,
        copy(smld.metadata)
    )
end

# The get_dimers function has been moved to analysis.jl