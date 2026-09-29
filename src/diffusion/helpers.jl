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
    dimerize(e1::DiffusingEmitter2D, e2::DiffusingEmitter2D, d_dimer::Float64; anchor::Union{Nothing,Int}=nothing)

Create two new emitters in dimer state from two monomers.

# Arguments
- `e1::DiffusingEmitter2D`: First emitter
- `e2::DiffusingEmitter2D`: Second emitter
- `d_dimer::Float64`: Dimer separation distance in microns
- `anchor::Union{Nothing,Int}=nothing`: `track_id` of the emitter that keeps its position; the other is placed
  `d_dimer` from it along the axis between them. `nothing` snaps both to the midpoint ± `d_dimer/2`.

# Returns
- `Tuple{DiffusingEmitter2D, DiffusingEmitter2D}`: Two new emitters in dimer state
"""
function dimerize(e1::DiffusingEmitter2D{T}, e2::DiffusingEmitter2D{T}, d_dimer::Float64; anchor::Union{Nothing,Int}=nothing) where T <: AbstractFloat
    # Calculate center of mass
    com_x = (e1.x + e2.x) / 2
    com_y = (e1.y + e2.y) / 2
    
    # Calculate orientation
    ϕ = angle(e1, e2)
    r = d_dimer / 2
    
    # Calculate new positions
    dx = r * cos(ϕ)
    dy = r * sin(ϕ)
    x1, y1, x2, y2 = com_x - dx, com_y - dy, com_x + dx, com_y + dy
    if anchor == e1.track_id
        x1, y1 = e1.x, e1.y
        x2, y2 = e1.x + 2dx, e1.y + 2dy
    elseif anchor == e2.track_id
        x2, y2 = e2.x, e2.y
        x1, y1 = e2.x - 2dx, e2.y - 2dy
    end
    
    # Create new dimer emitters
    d1 = DiffusingEmitter2D{T}(
        x1, y1,                  # Position
        e1.photons,              # Photons
        e1.timestamp,            # Timestamp
        e1.frame,                # Frame
        e1.dataset,              # Dataset
        e1.track_id,              # ID
        :dimer,                  # State
        e2.track_id                    # Partner ID
    )
    
    d2 = DiffusingEmitter2D{T}(
        x2, y2,                  # Position
        e2.photons,              # Photons
        e2.timestamp,            # Timestamp
        e2.frame,                # Frame
        e2.dataset,              # Dataset
        e2.track_id,              # ID
        :dimer,                  # State
        e1.track_id                    # Partner ID
    )
    
    return (d1, d2)
end

"""
    dimerize(e1::DiffusingEmitter3D, e2::DiffusingEmitter3D, d_dimer::Float64; anchor::Union{Nothing,Int}=nothing)

Create two new emitters in dimer state from two monomers in 3D.

# Arguments
- `e1::DiffusingEmitter3D`: First emitter
- `e2::DiffusingEmitter3D`: Second emitter
- `d_dimer::Float64`: Dimer separation distance in microns
- `anchor::Union{Nothing,Int}=nothing`: `track_id` of the emitter that keeps its position; the other is placed
  `d_dimer` from it along the axis between them. `nothing` snaps both to the midpoint ± `d_dimer/2`.

# Returns
- `Tuple{DiffusingEmitter3D, DiffusingEmitter3D}`: Two new emitters in dimer state
"""
function dimerize(e1::DiffusingEmitter3D{T}, e2::DiffusingEmitter3D{T}, d_dimer::Float64; anchor::Union{Nothing,Int}=nothing) where T <: AbstractFloat
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
    x1, y1, z1 = com_x - dx, com_y - dy, com_z - dz
    x2, y2, z2 = com_x + dx, com_y + dy, com_z + dz
    if anchor == e1.track_id
        x1, y1, z1 = e1.x, e1.y, e1.z
        x2, y2, z2 = e1.x + 2dx, e1.y + 2dy, e1.z + 2dz
    elseif anchor == e2.track_id
        x2, y2, z2 = e2.x, e2.y, e2.z
        x1, y1, z1 = e2.x - 2dx, e2.y - 2dy, e2.z - 2dz
    end
    
    # Create new dimer emitters
    d1 = DiffusingEmitter3D{T}(
        x1, y1, z1,                          # Position
        e1.photons,                          # Photons
        e1.timestamp,                        # Timestamp
        e1.frame,                            # Frame
        e1.dataset,                          # Dataset
        e1.track_id,                         # ID
        :dimer,                              # State
        e2.track_id                          # Partner ID
    )
    
    d2 = DiffusingEmitter3D{T}(
        x2, y2, z2,                          # Position
        e2.photons,                          # Photons
        e2.timestamp,                        # Timestamp
        e2.frame,                            # Frame
        e2.dataset,                          # Dataset
        e2.track_id,                         # ID
        :dimer,                              # State
        e1.track_id                          # Partner ID
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
    
    partner_idx = findfirst(em -> em.track_id == e.partner_id, emitters)
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
        e.track_id,          # ID
        :monomer,           # State
        nothing             # Partner ID
    )
    
    m2 = DiffusingEmitter2D{T}(
        partner.x, partner.y,  # Position
        partner.photons,       # Photons
        partner.timestamp,     # Timestamp
        partner.frame,         # Frame
        partner.dataset,       # Dataset
        partner.track_id,       # ID
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
    
    partner_idx = findfirst(em -> em.track_id == e.partner_id, emitters)
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
        e.track_id,         # ID
        :monomer,           # State
        nothing             # Partner ID
    )
    
    m2 = DiffusingEmitter3D{T}(
        partner.x, partner.y, partner.z,  # Position
        partner.photons,                  # Photons
        partner.timestamp,                # Timestamp
        partner.frame,                    # Frame
        partner.dataset,                  # Dataset
        partner.track_id,                 # ID
        :monomer,                         # State
        nothing                           # Partner ID
    )
    
    return (m1, m2)
end

# Relative margin on the separation of freshly dissociated partners, so that they end up
# strictly beyond `r_react` after rounding (also used by the DiffusionSMLMConfig box check).
const UNBIND_MARGIN = 1e-9

"""
    _anchor(params, e1, e2, D1, D2) -> Union{Nothing,Int}

The `track_id` of the pair member that keeps its position, or `nothing` when the pair is not
anchored. A pair is anchored only when `params.pair_mobility == :min` and one member is immobile
(`D == 0`); that member is the anchor, and if both are immobile the one with the smaller `track_id`.
"""
function _anchor(params, e1, e2, D1::Real, D2::Real)
    (params.pair_mobility == :min && (D1 == 0 || D2 == 0)) || return nothing
    return D1 == 0 && (D2 != 0 || e1.track_id < e2.track_id) ? e1.track_id : e2.track_id
end

# New positions (tuples) for the two members of `_unbind`; `p1`, `p2` are their coordinates.
function _unbind_positions(p1::NTuple{N,Float64}, p2::NTuple{N,Float64}, params, anchor, id1) where N
    s = params.r_react * (1 + UNBIND_MARGIN)
    d = p2 .- p1
    n = sqrt(sum(abs2, d))
    u = n > 0 ? d ./ n : ntuple(k -> k == 1 ? 1.0 : 0.0, N)
    reflecting = params.boundary == "reflecting"
    box = params.box_size
    if anchor !== nothing
        # Flip the axes on which the placed member would leave a reflecting box
        function place(a, sign)
            q = a .+ (sign * s) .* u
            reflecting || return q
            return ntuple(k -> 0 <= q[k] <= box ? q[k] : a[k] - (sign * s) * u[k], N)
        end
        return anchor == id1 ? (p1, place(p1, 1)) : (place(p2, -1), p2)
    end
    c = (p1 .+ p2) ./ 2
    if reflecting
        c = ntuple(k -> clamp(c[k], (s / 2) * abs(u[k]), box - (s / 2) * abs(u[k])), N)
    end
    return (c .- (s / 2) .* u, c .+ (s / 2) .* u)
end

"""
    _unbind(m1, m2, params, D1, D2)

Place two freshly dissociated monomers at least `params.r_react` apart along their pair axis,
so that the formation check (`can_dimerize`) cannot re-capture them at the next step only
because `d_dimer < r_react`. Partners already `r_react` or more apart are returned unchanged.
Under `pair_mobility = :min` an immobile member (both immobile: the lower `track_id`, see
`_anchor`) keeps its position; otherwise the pair midpoint is kept. In a reflecting box a member
that would leave the box is mirrored across the anchor on each axis it would leave, and a
symmetric pair's midpoint is shifted inward just far enough that both members are inside, so
`apply_boundary` cannot fold a member back within `r_react`; this needs `box_size > 2 r_react`.
Draws no random numbers and changes positions only.

# Arguments
- `m1, m2`: The monomers returned by `dissociate`
- `params`: Simulation parameters (`r_react`, `pair_mobility`, `box_size`, `boundary`)
- `D1, D2::Float64`: Monomer diffusion coefficients of `m1` and `m2`

# Returns
- `Tuple`: The two monomers, moved apart if they were closer than `r_react`
"""
function _unbind(m1::DiffusingEmitter2D{T}, m2::DiffusingEmitter2D{T}, params, D1::Real, D2::Real) where T <: AbstractFloat
    distance(m1, m2) >= params.r_react && return (m1, m2)
    anchor = _anchor(params, m1, m2, D1, D2)
    q1, q2 = _unbind_positions((Float64(m1.x), Float64(m1.y)), (Float64(m2.x), Float64(m2.y)), params, anchor, m1.track_id)
    mk(e, q) = DiffusingEmitter2D{T}(q[1], q[2], e.photons, e.timestamp, e.frame, e.dataset, e.track_id, e.state, e.partner_id)
    anchor == m1.track_id && return (m1, mk(m2, q2))
    anchor == m2.track_id && return (mk(m1, q1), m2)
    return (mk(m1, q1), mk(m2, q2))
end

function _unbind(m1::DiffusingEmitter3D{T}, m2::DiffusingEmitter3D{T}, params, D1::Real, D2::Real) where T <: AbstractFloat
    distance(m1, m2) >= params.r_react && return (m1, m2)
    anchor = _anchor(params, m1, m2, D1, D2)
    q1, q2 = _unbind_positions((Float64(m1.x), Float64(m1.y), Float64(m1.z)), (Float64(m2.x), Float64(m2.y), Float64(m2.z)), params, anchor, m1.track_id)
    mk(e, q) = DiffusingEmitter3D{T}(q[1], q[2], q[3], e.photons, e.timestamp, e.frame, e.dataset, e.track_id, e.state, e.partner_id)
    anchor == m1.track_id && return (m1, mk(m2, q2))
    anchor == m2.track_id && return (mk(m1, q1), m2)
    return (mk(m1, q1), mk(m2, q2))
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
        e.track_id,         # ID
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
        e.track_id,           # ID
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
        e1.track_id,             # ID
        e1.state,                # State
        e1.partner_id            # Partner ID
    )
    
    d2 = DiffusingEmitter2D{T}(
        com_x + dx, com_y + dy,  # Updated position
        e2.photons,              # Photons
        e2.timestamp + dt,       # Updated timestamp
        e2.frame,                # Frame
        e2.dataset,              # Dataset
        e2.track_id,             # ID
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
        e1.track_id,                         # ID
        e1.state,                            # State
        e1.partner_id                        # Partner ID
    )
    
    d2 = DiffusingEmitter3D{T}(
        com_x + dx, com_y + dy, com_z + dz,  # Updated position
        e2.photons,                          # Photons
        e2.timestamp + dt,                   # Updated timestamp
        e2.frame,                            # Frame
        e2.dataset,                          # Dataset
        e2.track_id,                         # ID
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
            e.track_id,        # ID
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
            e.track_id,           # ID
            e.state,              # State
            e.partner_id          # Partner ID
        )
    else
        return e
    end
end

# SMLD conversion utilities
"""
    create_smld(emitters::Vector{<:AbstractDiffusingEmitter}, camera::AbstractCamera, params::DiffusionSMLMConfig)

Convert a collection of diffusing emitters to a BasicSMLD object.

# Arguments
- `emitters::Vector{<:AbstractDiffusingEmitter}`: Collection of emitters from simulation
- `camera::AbstractCamera`: Camera model for imaging
- `params::DiffusionSMLMConfig`: Simulation parameters
- `track_D::Dict{Int,Float64}`: Per-track monomer diffusion coefficients, stored as `metadata["monomer_D"]`
- `track_class::Dict{Int,Int}`: Per-track mobility class (index into `monomer_mobility`), stored as `metadata["monomer_class"]`
- `γ::Union{Nothing,Real}=nothing`: Emission rate, photons/s, stored as `metadata["γ"]`; no key when `nothing`
- `n_frames::Union{Nothing,Int}=nothing`: Frame count of the movie; the largest frame present when `nothing`

# Returns
- `BasicSMLD`: SMLD containing all emitters for further analysis or visualization
"""
function create_smld(emitters::Vector{<:AbstractDiffusingEmitter}, camera::AbstractCamera, params::DiffusionSMLMConfig;
                     track_D::Dict{Int,Float64}=Dict{Int,Float64}(),
                     track_class::Dict{Int,Int}=Dict{Int,Int}(),
                     γ::Union{Nothing,Real}=nothing, n_frames::Union{Nothing,Int}=nothing)
    # Determine max frame number
    max_frame = n_frames !== nothing ? n_frames : isempty(emitters) ? 0 : maximum(e -> e.frame, emitters)
    
    # Create metadata
    metadata = Dict{String,Any}(
        "simulation_type" => "diffusion",
        "simulation_parameters" => params,
        "camera_framerate" => params.camera_framerate,
        "camera_exposure" => params.camera_exposure,
        "n_substeps" => substeps_per_frame(params)[1],
        "monomer_D" => track_D,
        "monomer_class" => track_class,
        "pair_mobility" => params.pair_mobility
    )
    γ === nothing || (metadata["γ"] = γ)
    
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