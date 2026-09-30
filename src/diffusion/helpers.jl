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

# Separation of freshly dissociated partners with coordinates of type T: r_react plus a margin that
# survives rounding each coordinate to T (|coordinate| ≤ box_size + r_react), at least a relative 1e-9
# (for Float64 the relative term dominates unless box_size exceeds about 5e5 × r_react). Internal.
_unbind_spacing(params, ::Type{T}) where {T<:AbstractFloat} =
    max(params.r_react * (1 + 1e-9),
        params.r_react + 8 * Float64(eps(T)) * max(1.0, params.box_size + params.r_react))

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

# New positions of the two members of `_unbind` in the coordinate type T, by the placement rule
# (dev/outputs/placement-rule.md, #39): `(ok, q1, q2)`, `ok` false when the construction cannot fit.
# `p1`, `p2` are their coordinates and `s` their separation (`_unbind_spacing`); under periodic
# boundaries the other member is taken at its minimum image. The anchor's own entry is its position.
function _unbind_positions(p1::NTuple{N,Float64}, p2::NTuple{N,Float64}, params, anchor, id1, s::Float64,
                           ::Type{T}) where {N,T<:AbstractFloat}
    box, reflecting = params.box_size, params.boundary == "reflecting"
    near(p, r) = reflecting ? p : ntuple(k -> p[k] - box * round((p[k] - r[k]) / box), Val(N))
    if anchor !== nothing
        a, m = anchor == id1 ? (p1, p2) : (p2, p1)
        ok, q = _place_from(a, _axis(a, near(m, a)), s, box, reflecting, T)
        at = ntuple(k -> T(a[k]), Val(N))
        return anchor == id1 ? (ok, at, q) : (ok, q, at)
    end
    p2 = near(p2, p1)
    u, c = _axis(p1, p2), (p1 .+ p2) ./ 2
    reflecting && return _place_centered(c, u, s, box, T)
    ok = all(ntuple(k -> s * abs(u[k]) <= box / 2, Val(N)))
    return ok, ntuple(k -> _inbox(T, mod(c[k] - (s / 2) * u[k], box), box), Val(N)),
           ntuple(k -> _inbox(T, mod(c[k] + (s / 2) * u[k], box), box), Val(N))
end

"""
    _unbind(m1, m2, params, D1, D2)

Place two freshly dissociated monomers at least `params.r_react` apart along their pair axis (the
minimum-image distance under periodic boundaries), so that the formation check (`can_dimerize`) cannot
re-capture them at the next step only because `d_dimer < r_react`. Partners already `r_react` or more
apart are returned unchanged. Under `pair_mobility = :min` an immobile member (both immobile: the lower
`track_id`, see `_anchor`) keeps its position and the other is placed from it (mirrored across it on each
axis it would leave a reflecting box); otherwise, including every pair under `:fixed`, both move apart
about their midpoint, shifted inward just far enough that both are inside a reflecting box, or wrapped
under periodic boundaries. Every placed coordinate is inside the box in the coordinate type, and the
separation carries a margin in that precision (`_unbind_spacing`), so Float32 partners are also beyond
`r_react` after conversion. When that placement cannot fit along the pair's axis (see
dev/outputs/placement-rule.md), the partners are left where `dissociate` put them (0.7.2's behaviour),
with a one-time warning. Draws no random numbers and changes positions only.

# Arguments
- `m1, m2`: The monomers returned by `dissociate`
- `params`: Simulation parameters (`r_react`, `pair_mobility`, `box_size`, `boundary`)
- `D1, D2::Float64`: Monomer diffusion coefficients of `m1` and `m2`

# Returns
- `Tuple`: The two monomers, moved apart if they were closer than `r_react`
"""
function _unbind(m1::E, m2::E, params, D1::Real, D2::Real) where {E<:AbstractDiffusingEmitter}
    distance(m1, _near(m2, m1, params)) >= params.r_react && return (m1, m2)
    T = typeof(m1.x)
    anchor = _anchor(params, m1, m2, D1, D2)
    ok, q1, q2 = _unbind_positions(_pos(m1), _pos(m2), params, anchor, m1.track_id, _unbind_spacing(params, T), T)
    if !ok
        @warn "box_size=$(params.box_size) is too small to place dissociated partners r_react=$(params.r_react) apart along their axis; they stay where they are (0.7.2's behaviour) and may re-form at once" maxlog=1
        return (m1, m2)
    end
    anchor == m1.track_id && return (m1, _at(m2, q2))
    anchor == m2.track_id && return (_at(m1, q1), m2)
    return (_at(m1, q1), _at(m2, q2))
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

# Rigid-pair placement (dev/outputs/placement-rule.md), shared by pair formation and bound motion.
# Positions are computed in Float64 from the partners' coordinates and returned in the coordinate
# type T, inside [0, _top(T, box)] on every axis. Internal.

# The largest value of T that is <= box, and x converted to T and clamped into [0, that value]
_top(::Type{T}, box::Float64) where {T<:AbstractFloat} = T(box) <= box ? T(box) : prevfloat(T(box))
_inbox(::Type{T}, x::Float64, box::Float64) where {T<:AbstractFloat} = clamp(T(x), zero(T), _top(T, box))

_coords(e::DiffusingEmitter2D) = (e.x, e.y)
_coords(e::DiffusingEmitter3D) = (e.x, e.y, e.z)
_pos(e::AbstractDiffusingEmitter) = Float64.(_coords(e))
_at(e::DiffusingEmitter2D{T}, p::NTuple{2,T}) where {T} =
    DiffusingEmitter2D{T}(p[1], p[2], e.photons, e.timestamp, e.frame, e.dataset, e.track_id, e.state, e.partner_id)
_at(e::DiffusingEmitter3D{T}, p::NTuple{3,T}) where {T} =
    DiffusingEmitter3D{T}(p[1], p[2], p[3], e.photons, e.timestamp, e.frame, e.dataset, e.track_id, e.state,
                          e.partner_id)
_inside(e::AbstractDiffusingEmitter, box::Float64) = all(c -> 0 <= c <= box, _coords(e))

# `e` at the minimum image of its offset from `ref` under periodic boundaries, so a pair's center and
# axis are those of its bond, not of its wrapped coordinates; unchanged under reflecting boundaries
function _near(e::AbstractDiffusingEmitter, ref::AbstractDiffusingEmitter, params)
    params.boundary == "periodic" || return e
    L = params.box_size
    return _at(e, map((c, r) -> typeof(c)(c - L * round((c - r) / L)), _coords(e), _coords(ref)))
end

# x folded into [lo, hi] as often as it crosses either end (a triangle wave); lo when hi <= lo
function _fold(x::Float64, lo::Float64, hi::Float64)
    hi <= lo && return lo
    w = hi - lo
    y = mod(x - lo, 2w)
    return lo + (y <= w ? y : 2w - y)
end

# Unit vector from p to q; (1, 0, ...) when they coincide
function _axis(p::NTuple{N,Float64}, q::NTuple{N,Float64}) where {N}
    v = q .- p
    n = sqrt(sum(abs2, v))
    return n > 0 ? v ./ n : ntuple(k -> k == 1 ? 1.0 : 0.0, Val(N))
end

# The point `s` from the fixed point `a` along `u`. Reflecting: each axis that would leave the box
# [0, box] is mirrored across `a`, and `ok` is false when some axis fits neither way. Periodic: wrapped,
# and `ok` is false unless every s|u_k| <= box/2, so the offset is its own minimum image.
function _place_from(a::NTuple{N,Float64}, u::NTuple{N,Float64}, s::Float64, box::Float64, reflecting::Bool,
                     ::Type{T}) where {N,T<:AbstractFloat}
    fits(x) = 0 <= x <= box
    ok = all(ntuple(k -> reflecting ? fits(a[k] + s * u[k]) || fits(a[k] - s * u[k]) : s * abs(u[k]) <= box / 2, Val(N)))
    q = ntuple(Val(N)) do k
        x = a[k] + s * u[k]
        _inbox(T, reflecting ? (fits(x) ? x : a[k] - s * u[k]) : mod(x, box), box)
    end
    return ok, q
end

# Two points `s` apart along `u`, centered on `c`, in a reflecting box [0, box]: the center keeps
# (s/2)|u_k| from each wall, clamped there when placing, or folded into that interval as often as it
# crosses when moving (`reflect = true`); `ok` is false when some s|u_k| exceeds the box.
function _place_centered(c::NTuple{N,Float64}, u::NTuple{N,Float64}, s::Float64, box::Float64,
                         ::Type{T}; reflect::Bool=false) where {N,T<:AbstractFloat}
    h = ntuple(k -> (s / 2) * abs(u[k]), Val(N))
    ok = all(ntuple(k -> 2h[k] <= box, Val(N)))
    m = ntuple(k -> reflect ? _fold(c[k], h[k], box - h[k]) : clamp(c[k], h[k], box - h[k]), Val(N))
    return ok, ntuple(k -> _inbox(T, m[k] - (s / 2) * u[k], box), Val(N)), ntuple(k -> _inbox(T, m[k] + (s / 2) * u[k], box), Val(N))
end

# A pair just formed from monomers `e1`, `e2` (`d1`, `d2` from `dimerize`), by the placement rule: an
# anchored pair keeps the anchor and places the other partner `d_dimer` from it, or leaves it where it
# was when that cannot fit; a mobile pair in a reflecting box with an end outside moves to its midpoint
# shifted inward just enough (each end reflected on its own when it cannot fit); a mobile pair under
# periodic boundaries is left as in 0.7.1. Internal.
function _place_pair(d1::E, d2::E, e1::E, e2::E, anchor::Union{Nothing,Int},
                     params::DiffusionSMLMConfig) where {E<:AbstractDiffusingEmitter}
    T = typeof(d1.x)
    box, reflecting = params.box_size, params.boundary == "reflecting"
    if anchor !== nothing
        fixed, mover = anchor == e1.track_id ? (e1, e2) : (e2, e1)
        a = _pos(fixed)
        ok, q = _place_from(a, _axis(a, _pos(mover)), params.d_dimer, box, reflecting, T)
        p = ok ? q : _coords(mover)
        return anchor == e1.track_id ? (d1, _at(d2, p)) : (_at(d1, p), d2)
    end
    (!reflecting || (_inside(d1, box) && _inside(d2, box))) && return d1, d2
    p1, p2 = _pos(e1), _pos(e2)
    ok, q1, q2 = _place_centered((p1 .+ p2) ./ 2, _axis(p1, p2), params.d_dimer, box, T)
    ok && return _at(d1, q1), _at(d2, q2)
    return apply_boundary(d1, box, params.boundary), apply_boundary(d2, box, params.boundary)
end

# A mobile pair after a bound step (`d1`, `d2` from `diffuse_dimer` of the partner's minimum image,
# `_near`): in a reflecting box with an end outside, its center folds off the walls moved in by each
# end's half-extent, keeping orientation and bond length (each end reflected on its own when it cannot
# fit); under periodic boundaries each end is wrapped. Internal.
function _move_pair(d1::E, d2::E, params::DiffusionSMLMConfig) where {E<:AbstractDiffusingEmitter}
    box = params.box_size
    params.boundary == "reflecting" || return apply_boundary(d1, box, params.boundary), apply_boundary(d2, box, params.boundary)
    _inside(d1, box) && _inside(d2, box) && return d1, d2
    p1, p2 = _pos(d1), _pos(d2)
    ok, q1, q2 = _place_centered((p1 .+ p2) ./ 2, _axis(p1, p2), params.d_dimer, box, typeof(d1.x); reflect=true)
    ok && return _at(d1, q1), _at(d2, q2)
    return apply_boundary(d1, box, params.boundary), apply_boundary(d2, box, params.boundary)
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
    
    # Inside [0, box_size] in T: a Float32 coordinate can round past the wall
    new_x, new_y = _inbox(T, Float64(new_x), box_size), _inbox(T, Float64(new_y), box_size)

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
    
    # Inside [0, box_size] in T: a Float32 coordinate can round past the wall
    new_x, new_y, new_z = _inbox(T, Float64(new_x), box_size), _inbox(T, Float64(new_y), box_size),
                          _inbox(T, Float64(new_z), box_size)

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
- `params::DiffusionSMLMConfig`: Simulation parameters, stored as `metadata["simulation_parameters"]`; its `dt`
  is also stored as `metadata["dt"]`, a snapshot that later edits of the mutable config do not change
- `track_D::Dict{Int,Float64}`: Per-track monomer diffusion coefficients, stored as `metadata["monomer_D"]`; the
  config's mixture is stored as `metadata["monomer_mobility"]`, a copy that later edits of the config do not change
- `track_class::Dict{Int,Int}`: Per-track mobility class (index into `monomer_mobility`), stored as `metadata["monomer_class"]`
- `γ::Union{Nothing,Real}=nothing`: Emission rate, photons/s, stored as `metadata["γ"]`; no key when `nothing`
- `rate_source::Union{Nothing,String}=nothing`: How the rate was set (`"γ"`, `"photons"` or `"default"`), stored as
  `metadata["rate_source"]`; no key when `nothing`
- `n_frames::Union{Nothing,Int}=nothing`: Frame count of the movie; the largest frame present when `nothing`

# Returns
- `BasicSMLD`: SMLD containing all emitters for further analysis or visualization
"""
function create_smld(emitters::Vector{<:AbstractDiffusingEmitter}, camera::AbstractCamera, params::DiffusionSMLMConfig;
                     track_D::Dict{Int,Float64}=Dict{Int,Float64}(),
                     track_class::Dict{Int,Int}=Dict{Int,Int}(),
                     γ::Union{Nothing,Real}=nothing, rate_source::Union{Nothing,String}=nothing,
                     n_frames::Union{Nothing,Int}=nothing)
    # Determine max frame number
    max_frame = n_frames !== nothing ? n_frames : isempty(emitters) ? 0 : maximum(e -> e.frame, emitters)
    
    # Create metadata
    metadata = Dict{String,Any}(
        "simulation_type" => "diffusion",
        "simulation_parameters" => params,
        "dt" => params.dt,
        "camera_framerate" => params.camera_framerate,
        "camera_exposure" => params.camera_exposure,
        "n_substeps" => substeps_per_frame(params)[1],
        "monomer_D" => track_D,
        "monomer_class" => track_class,
        "monomer_mobility" => copy(params.monomer_mobility),
        "pair_mobility" => params.pair_mobility
    )
    γ === nothing || (metadata["γ"] = γ)
    rate_source === nothing || (metadata["rate_source"] = rate_source)
    
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