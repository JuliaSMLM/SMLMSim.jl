"""
    DiffusionSMLMConfig <: SMLMSimParams

Parameters for diffusion-based SMLM simulation using Smoluchowski dynamics.

# Fields
- `density::Float64`: number density (molecules/μm²)
- `box_size::Float64`: simulation box size (μm)
- `diff_monomer::Float64`: monomer diffusion coefficient (μm²/s)
- `diff_dimer::Float64`: dimer diffusion coefficient (μm²/s)
- `diff_dimer_rot::Float64`: dimer rotational diffusion coefficient (rad²/s)
- `k_off::Float64`: dimer dissociation rate (s⁻¹)
- `r_react::Float64`: reaction radius (μm)
- `d_dimer::Float64`: monomer separation in dimer (μm)
- `dt::Float64`: physics step (s); also sets the sub-steps per frame (motion blur):
  camera_exposure and 1/camera_framerate should be integer multiples of dt (otherwise
  `simulate` rounds to the nearest step count and warns)
- `t_max::Float64`: total simulation time (s); rounded down to whole frame periods
  (`t_max` shorter than one frame period simulates one frame, with a warning)
- `ndims::Int`: number of dimensions (2 or 3)
- `boundary::String`: boundary condition type ("periodic" or "reflecting")
- `camera_framerate::Float64`: camera frames per second (Hz)
- `camera_exposure::Float64`: camera exposure time per frame (s)
- `monomer_mobility::Vector{Tuple{Float64,Float64}}`: optional mixture of monomer mobility
  populations, each entry `(fraction, D)` with fractions summing to 1 and `D` in μm²/s.
  Each molecule draws its population once, at initialization (or when it enters through
  `starting_conditions`), and keeps its own `D` as a monomer, including after dissociation.
  Empty (default) means every monomer uses `diff_monomer`. Dimers always use `diff_dimer`.
  The drawn values are stored in `smld.metadata["monomer_D"]` (track_id => D).

Photons: `simulate` takes `γ`, the emission rate in photons/s. Each of the
`n_sub` records of a frame (`substeps_per_frame`) carries `γ·dt`, so a frame holds `γ·n_sub·dt`
photons: `γ·camera_exposure` when `camera_exposure` is an integer multiple of `dt` and not longer
than the frame period `1/camera_framerate`. There is no field for γ; it is a keyword of `simulate`.

# Examples
```julia
# Default parameters
params = DiffusionSMLMConfig()

# Custom parameters
params = DiffusionSMLMConfig(
    density = 1.0,           # 1 molecule per μm²
    box_size = 20.0,         # 20μm × 20μm box
    diff_monomer = 0.2,      # 0.2 μm²/s
    diff_dimer = 0.1,        # 0.1 μm²/s
    diff_dimer_rot = 0.8,    # 0.8 rad²/s
    k_off = 0.1,             # 0.1 s⁻¹
    r_react = 0.02,          # 20nm reaction radius
    d_dimer = 0.06,          # 60nm dimer separation
    dt = 0.005,              # 5ms time step
    t_max = 20.0,            # 20s simulation
    ndims = 3,               # 3D simulation
    boundary = "reflecting", # reflecting boundaries
    camera_framerate = 20.0, # 20 frames per second
    camera_exposure = 0.04   # 40ms exposure per frame
)

# Mixed mobility: 85% immobile, 5% slow, 10% mobile monomers
params = DiffusionSMLMConfig(
    monomer_mobility = [(0.85, 0.0), (0.05, 0.08), (0.10, 0.38)]
)
```
"""
Base.@kwdef mutable struct DiffusionSMLMConfig <: SMLMSimParams
    density::Float64 = 1.0
    box_size::Float64 = 10.0
    diff_monomer::Float64 = 0.1
    diff_dimer::Float64 = 0.05
    diff_dimer_rot::Float64 = 0.5
    k_off::Float64 = 0.2
    r_react::Float64 = 0.01
    d_dimer::Float64 = 0.05
    dt::Float64 = 0.01
    t_max::Float64 = 10.0
    ndims::Int = 2
    boundary::String = "periodic"
    
    # Camera imaging parameters
    camera_framerate::Float64 = 10.0     # Frames per second
    camera_exposure::Float64 = 0.1       # Exposure time in seconds

    # Monomer mobility mixture: (fraction, D) entries; empty means diff_monomer for all
    monomer_mobility::Vector{Tuple{Float64,Float64}} = Tuple{Float64,Float64}[]
    
    function DiffusionSMLMConfig(
        density, box_size, diff_monomer, diff_dimer, diff_dimer_rot,
        k_off, r_react, d_dimer, dt, t_max, ndims, boundary,
        camera_framerate, camera_exposure, monomer_mobility
    )
        # Input validation
        if density <= 0
            throw(ArgumentError("Density must be positive"))
        end
        if box_size <= 0
            throw(ArgumentError("Box size must be positive"))
        end
        if diff_monomer < 0
            throw(ArgumentError("Monomer diffusion coefficient must be non-negative"))
        end
        if diff_dimer < 0
            throw(ArgumentError("Dimer diffusion coefficient must be non-negative"))
        end
        if diff_dimer_rot < 0
            throw(ArgumentError("Dimer rotational diffusion coefficient must be non-negative"))
        end
        if k_off < 0
            throw(ArgumentError("Dissociation rate must be non-negative"))
        end
        if r_react <= 0
            throw(ArgumentError("Reaction radius must be positive"))
        end
        if d_dimer <= 0
            throw(ArgumentError("Dimer separation must be positive"))
        end
        if dt <= 0
            throw(ArgumentError("Time step must be positive"))
        end
        if t_max <= 0
            throw(ArgumentError("Maximum simulation time must be positive"))
        end
        if ndims != 2 && ndims != 3
            throw(ArgumentError("Number of dimensions must be 2 or 3"))
        end
        if boundary != "periodic" && boundary != "reflecting"
            throw(ArgumentError("Boundary condition must be 'periodic' or 'reflecting'"))
        end
        # Camera parameter validation
        if camera_framerate <= 0
            throw(ArgumentError("Camera framerate must be positive"))
        end
        if camera_exposure <= 0
            throw(ArgumentError("Camera exposure time must be positive"))
        end
        # Monomer mobility mixture validation
        if !isempty(monomer_mobility)
            if any(m -> !(m[1] > 0), monomer_mobility)
                throw(ArgumentError("Monomer mobility fractions must be positive"))
            end
            if any(m -> !(isfinite(m[2]) && m[2] >= 0), monomer_mobility)
                throw(ArgumentError("Monomer mobility diffusion coefficients must be non-negative"))
            end
            if !(abs(sum(m -> m[1], monomer_mobility) - 1) <= 1e-9)
                throw(ArgumentError("Monomer mobility fractions must sum to 1"))
            end
        end
        
        new(density, box_size, diff_monomer, diff_dimer, diff_dimer_rot,
            k_off, r_react, d_dimer, dt, t_max, ndims, boundary,
            camera_framerate, camera_exposure, monomer_mobility)
    end
end

# Positional construction without a mobility mixture (the 0.7 signature)
DiffusionSMLMConfig(density, box_size, diff_monomer, diff_dimer, diff_dimer_rot,
                    k_off, r_react, d_dimer, dt, t_max, ndims, boundary,
                    camera_framerate, camera_exposure) =
    DiffusionSMLMConfig(density, box_size, diff_monomer, diff_dimer, diff_dimer_rot,
                        k_off, r_react, d_dimer, dt, t_max, ndims, boundary,
                        camera_framerate, camera_exposure, Tuple{Float64,Float64}[])

"""
    substeps_per_frame(params::DiffusionSMLMConfig) -> (n_sub, steps_per_frame)

Number of physics steps `dt` inside one exposure (`n_sub`) and in one frame period
(`steps_per_frame`). Internal. Timing that is not an integer multiple of `dt` is rounded
to the nearest step count, and an exposure longer than the frame period is capped at the
frame period; each case warns (once) with the rule applied.
"""
function substeps_per_frame(params::DiffusionSMLMConfig)
    r_sub = params.camera_exposure / params.dt
    r_frame = 1 / (params.camera_framerate * params.dt)
    n_sub = max(1, round(Int, r_sub))
    steps_per_frame = max(1, round(Int, r_frame))
    if abs(r_sub - n_sub) > 1e-9 * max(1, r_sub)
        @warn "camera_exposure=$(params.camera_exposure) is not an integer multiple of dt=$(params.dt); using n_sub=$n_sub sub-steps per frame (effective exposure $(n_sub * params.dt) s)" maxlog=1
    end
    if abs(r_frame - steps_per_frame) > 1e-9 * max(1, r_frame)
        @warn "the frame period 1/camera_framerate=$(1 / params.camera_framerate) is not an integer multiple of dt=$(params.dt); using $steps_per_frame steps per frame (effective frame period $(steps_per_frame * params.dt) s)" maxlog=1
    end
    if n_sub > steps_per_frame
        @warn "camera_exposure=$(params.camera_exposure) is longer than the frame period 1/camera_framerate=$(1 / params.camera_framerate); capping the exposure at the frame period ($steps_per_frame sub-steps); use camera_exposure ≤ 1/camera_framerate" maxlog=1
        n_sub = steps_per_frame
    end
    return n_sub, steps_per_frame
end

"""
    draw_monomer_D(params::DiffusionSMLMConfig) -> Float64

Draw a monomer diffusion coefficient from `params.monomer_mobility`
(`params.diff_monomer` when the mixture is empty). Internal.
"""
function draw_monomer_D(params::DiffusionSMLMConfig)
    isempty(params.monomer_mobility) && return params.diff_monomer
    u = rand()
    c = 0.0
    for (frac, D) in params.monomer_mobility
        c += frac
        u < c && return D
    end
    return params.monomer_mobility[end][2]
end

"""
    initialize_emitters(params::DiffusionSMLMConfig, photons::Union{Nothing,Real}=nothing;
                        γ=nothing, override_count::Union{Nothing, Int}=nothing)

Create initial emitter positions for the simulation.

# Arguments
- `params::DiffusionSMLMConfig`: Simulation parameters
- `photons`: Deprecated (removed in 0.8.0). Photons per emitter record, as in 0.7.
  Pass `γ` instead.

# Keyword Arguments
- `γ`: Emission rate, photons/s; each emitter carries `γ·dt` photons (the photons of one
  sub-step), so a frame holds γ·camera_exposure photons. Passing both `photons` and `γ`
  throws `ArgumentError`. With neither, each emitter carries 1000.0 photons.
- `override_count::Union{Nothing, Int}=nothing`: Optional override for the number of molecules

# Returns
- `Vector{<:AbstractDiffusingEmitter}`: Vector of initialized emitters
"""
function initialize_emitters(params::DiffusionSMLMConfig, photons::Union{Nothing,Real}=nothing;
                             γ::Union{Nothing,Real}=nothing, override_count::Union{Nothing, Int}=nothing)
    photons !== nothing && γ !== nothing &&
        throw(ArgumentError("pass γ or the deprecated positional photons, not both"))
    if photons !== nothing
        Base.depwarn("the positional photons argument is deprecated and will be removed in 0.8.0; pass the keyword γ (photons/s)", :initialize_emitters; force=true)
        per_record = Float64(photons)
    elseif γ !== nothing
        per_record = Float64(γ) * params.dt
    else
        per_record = 1000.0
    end
    return build_emitters(params, per_record, override_count)
end

# Emitters at random positions, each carrying `photons` per record. Internal.
function build_emitters(params::DiffusionSMLMConfig, photons::Float64, override_count::Union{Nothing, Int})
    # Calculate number of molecules
    n_molecules = if override_count !== nothing
        override_count
    else
        round(Int, params.density * params.box_size^params.ndims)
    end
    
    # Create array for emitters
    if params.ndims == 2
        emitters = Vector{DiffusingEmitter2D{Float64}}(undef, n_molecules)
        
        # Create emitters at random positions
        for i in 1:n_molecules
            x = rand(Uniform(0, params.box_size))
            y = rand(Uniform(0, params.box_size))
            
            # Create emitter with initial properties
            emitters[i] = DiffusingEmitter2D{Float64}(
                x, y,                      # Position
                photons,                   # Photons per record
                0.0,                       # Initial timestamp
                1,                         # Initial frame
                1,                         # Dataset
                i,                         # track_id
                :monomer,                  # Initial state
                nothing                    # No partner initially
            )
        end
    else  # 3D
        emitters = Vector{DiffusingEmitter3D{Float64}}(undef, n_molecules)
        
        # Create emitters at random positions
        for i in 1:n_molecules
            x = rand(Uniform(0, params.box_size))
            y = rand(Uniform(0, params.box_size))
            z = rand(Uniform(0, params.box_size))
            
            # Create emitter with initial properties
            emitters[i] = DiffusingEmitter3D{Float64}(
                x, y, z,                   # Position
                photons,                   # Photons per record
                0.0,                       # Initial timestamp
                1,                         # Initial frame
                1,                         # Dataset
                i,                         # track_id
                :monomer,                  # Initial state
                nothing                    # No partner initially
            )
        end
    end
    
    return emitters
end

"""
    update_system(emitters::Vector{<:AbstractDiffusingEmitter}, params::DiffusionSMLMConfig, dt::Float64;
                  track_D=nothing)

Update all emitters based on Smoluchowski diffusion dynamics. Monomers diffuse with
`track_D[track_id]` when the track has an entry, otherwise with `params.diff_monomer`.

# Arguments
- `emitters::Vector{<:AbstractDiffusingEmitter}`: Current emitters state
- `params::DiffusionSMLMConfig`: Simulation parameters
- `dt::Float64`: Time step
- `track_D::Union{Nothing,Dict{Int,Float64}}=nothing`: Per-track monomer diffusion coefficients

# Returns
- `Vector{<:AbstractDiffusingEmitter}`: Updated emitters
"""
function update_system(emitters::Vector{<:AbstractDiffusingEmitter}, params::DiffusionSMLMConfig, dt::Float64;
                       track_D::Union{Nothing,Dict{Int,Float64}}=nothing)
    monomer_D(id) = track_D === nothing ? params.diff_monomer : get(track_D, id, params.diff_monomer)
    # Create new array for updated emitters
    new_emitters = Vector{eltype(emitters)}()
    
    # Create a set of IDs that have been processed
    processed = Set{Int}()
    
    # Process all emitters
    for (i, e1) in enumerate(emitters)
        e1.track_id in processed && continue
        
        if e1.state == :monomer
            # Check if this monomer forms a dimer with any other monomer
            found_dimer = false
            
            for (j, e2) in enumerate(emitters[i+1:end])
                e2.track_id in processed && continue
                e2.state == :monomer || continue
                
                if can_dimerize(e1, e2, params.r_react)
                    # Create new dimer pair
                    d1, d2 = dimerize(e1, e2, params.d_dimer)
                    push!(new_emitters, d1, d2)
                    push!(processed, e1.track_id, e2.track_id)
                    found_dimer = true
                    break
                end
            end
            
            # If didn't form a dimer, update as monomer
            if !found_dimer
                # Apply diffusion
                new_e = diffuse(e1, monomer_D(e1.track_id), dt)
                
                # Apply boundary conditions
                new_e = apply_boundary(new_e, params.box_size, params.boundary)
                
                push!(new_emitters, new_e)
                push!(processed, e1.track_id)
            end
        elseif e1.state == :dimer && !(e1.track_id in processed)
            # Check for dissociation
            if should_dissociate(e1, params.k_off, dt)
                # Find partner and create two new monomers
                m1, m2 = dissociate(e1, emitters)
                
                # Apply diffusion to each new monomer
                m1 = diffuse(m1, monomer_D(m1.track_id), dt)
                m2 = diffuse(m2, monomer_D(m2.track_id), dt)
                
                # Apply boundary conditions
                m1 = apply_boundary(m1, params.box_size, params.boundary)
                m2 = apply_boundary(m2, params.box_size, params.boundary)
                
                push!(new_emitters, m1, m2)
                push!(processed, e1.track_id, e1.partner_id)
            else
                # Find partner and update dimer
                partner_idx = findfirst(e -> e.track_id == e1.partner_id, emitters)
                if !isnothing(partner_idx) && !(emitters[partner_idx].track_id in processed)
                    e2 = emitters[partner_idx]
                    
                    # Apply dimer diffusion
                    d1, d2 = diffuse_dimer(
                        e1, e2, 
                        params.diff_dimer, 
                        params.diff_dimer_rot, 
                        params.d_dimer, 
                        dt
                    )
                    
                    # Apply boundary conditions
                    d1 = apply_boundary(d1, params.box_size, params.boundary)
                    d2 = apply_boundary(d2, params.box_size, params.boundary)
                    
                    push!(new_emitters, d1, d2)
                    push!(processed, e1.track_id, e2.track_id)
                else
                    @warn "dimer partner $(e1.partner_id) of track $(e1.track_id) not found; the emitter is dropped, as in 0.7" maxlog=1
                end
            end
        end
    end
    
    return new_emitters
end

"""
    restamp(e::AbstractDiffusingEmitter; photons, timestamp, frame, state, partner_id)

Copy of a diffusing emitter with the given fields replaced. Internal.
"""
restamp(e::DiffusingEmitter2D{T}; photons=e.photons, timestamp=e.timestamp, frame=e.frame,
        state=e.state, partner_id=e.partner_id) where T =
    DiffusingEmitter2D{T}(e.x, e.y, photons, timestamp, frame, e.dataset, e.track_id, state, partner_id)
restamp(e::DiffusingEmitter3D{T}; photons=e.photons, timestamp=e.timestamp, frame=e.frame,
        state=e.state, partner_id=e.partner_id) where T =
    DiffusingEmitter3D{T}(e.x, e.y, e.z, photons, timestamp, frame, e.dataset, e.track_id, state, partner_id)

"""
    add_camera_frame_emitters!(camera_emitters, emitters, time, frame_num, params)

Deprecated, removed in 0.8.0 (internal). Add emitters to camera frames when `time` falls
within the exposure window of frame `frame_num`; `simulate` no longer calls it.

# Arguments
- `camera_emitters::Vector{<:AbstractDiffusingEmitter}`: Collection of emitters for camera frames
- `emitters::Vector{<:AbstractDiffusingEmitter}`: Current emitters from simulation
- `time::Float64`: Current simulation time
- `frame_num::Int`: Current frame number
- `params::DiffusionSMLMConfig`: Simulation parameters

# Returns
- `Nothing`
"""
function add_camera_frame_emitters!(camera_emitters, emitters, time, frame_num, params::DiffusionSMLMConfig)
    Base.depwarn("add_camera_frame_emitters! is internal and deprecated; it will be removed in 0.8.0", :add_camera_frame_emitters!; force=true)
    # Check if this timepoint falls within a camera exposure window
    exposure_start = (frame_num - 1) / params.camera_framerate
    exposure_end = exposure_start + params.camera_exposure

    if time >= exposure_start && time <= exposure_end
        # Add all current emitters to the camera frame
        for e in emitters
            push!(camera_emitters, restamp(e; timestamp=time, frame=frame_num))
        end
    end

    return nothing
end

"""
    _record_frame!(camera_emitters, emitters, time, frame_num)

Record the current emitters as one sub-step record of camera frame `frame_num`.
Each record carries the live emitter's photons (`γ·dt` for a rate `γ`), so the `n_sub`
records of a frame sum to `γ·n_sub·dt`. `simulate` decides by integer step which
sub-steps to record. Internal.
"""
function _record_frame!(camera_emitters, emitters, time, frame_num)
    for e in emitters
        push!(camera_emitters, restamp(e; timestamp=time, frame=frame_num))
    end
    return nothing
end

"""
    simulate(params::DiffusionSMLMConfig;
             starting_conditions::Union{Nothing, SMLD, Vector{<:AbstractDiffusingEmitter}}=nothing,
             γ::Union{Nothing, Real}=nothing,
             photons::Union{Nothing, Real}=nothing,
             override_count::Union{Nothing, Int}=nothing,
             kwargs...)

Run a Smoluchowski diffusion simulation and return a BasicSMLD object
with emitters that have both frame number and timestamp information.

# Arguments
- `params::DiffusionSMLMConfig`: Simulation parameters

# Keyword Arguments
- `starting_conditions::Union{Nothing, SMLD, Vector{<:AbstractDiffusingEmitter}}=nothing`: Optional starting emitters.
  An SMLD resumes at its exact end state (`extract_end_state`) and carries each track's D and γ
  forward. A Vector keeps each emitter's own `photons`, gets fresh D draws, and is deduplicated
  to the latest record per track_id (with a warning) if track_ids repeat.
- `γ::Union{Nothing, Real}=nothing`: emission rate, photons/s (finite, ≥ 0); each of the
  n_sub records in a frame carries γ·dt, so a frame holds γ·n_sub·dt photons (γ·camera_exposure
  when the exposure is a whole number of steps and not capped at the frame period).
  An explicit γ restamps starting emitters to γ·dt. Default for new emitters: 1000 photons
  per record (γ = 1000/dt, 0.7's default; 0.8.0 will change the default to a fixed rate).
  `smld.metadata["rate_source"]` records how the rate was set: `"γ"`, `"photons"` or `"default"`.
  Without γ, an SMLD `starting_conditions` whose rate was set with γ continues at that rate (every
  molecule restamped to γ·dt at this `dt`); any other source (default, `photons`, a 0.7.1 SMLD, a
  Vector) keeps each molecule's photons per record, as in 0.7.1.
- `photons::Union{Nothing, Real}=nothing`: deprecated, removed in 0.8.0. Photons per record,
  as in 0.7; for new emitters `photons = p` is `γ = p/dt` with identical output. It is ignored with
  `starting_conditions`, as in 0.7, where `γ` restamps every molecule. Passing both `photons` and
  `γ` throws `ArgumentError`.
- `override_count::Union{Nothing, Int}=nothing`: Optional override for the number of molecules
- `camera::Union{Nothing, AbstractCamera}=nothing`: Camera model (default: IdealCamera with 100nm pixels)
  - If `nothing`, creates IdealCamera with dimensions matching box_size
  - Can specify SCMOSCamera for realistic noise modeling
- Any additional parameters are ignored (allows unified interface with other simulate methods)

# Warnings
Each warns once and then runs with the stated rule: `camera_exposure` or the frame period not an
integer multiple of `dt` (rounded to the nearest step count), `camera_exposure` longer than the
frame period (capped at the frame period), `t_max` under one frame period (one frame is
simulated), Vector `starting_conditions` with repeated track_id (latest record per track),
and dimers without a matching partner (converted to monomers).

# Returns
- `Tuple{BasicSMLD, SimInfo}`: (smld, info)
    - smld: SMLD object containing all emitters across all frames
    - info: SimInfo containing timing and simulation statistics

# Example
```julia
# Set up parameters with camera settings
params = DiffusionSMLMConfig(
    density = 0.5,           # molecules per μm²
    box_size = 10.0,         # μm
    camera_framerate = 20.0, # 20 fps
    camera_exposure = 0.04   # 40ms exposure
)

# Run basic simulation
smld, info = simulate(params)

# Run simulation with exactly 2 particles
smld, info = simulate(params; override_count=2)

# Use previous simulation state as starting conditions for a new simulation
smld_continued, info = simulate(params; starting_conditions=smld)
```
"""
function simulate(params::DiffusionSMLMConfig;
                 starting_conditions::Union{Nothing, SMLD, Vector{<:AbstractDiffusingEmitter}}=nothing,
                 γ::Union{Nothing, Real}=nothing,
                 photons::Union{Nothing, Real}=nothing,
                 override_count::Union{Nothing, Int}=nothing,
                 camera::Union{Nothing, AbstractCamera}=nothing,
                 kwargs...)

    start_time = time_ns()

    photons !== nothing && γ !== nothing &&
        throw(ArgumentError("pass γ or the deprecated photons, not both"))
    γ === nothing || (isfinite(γ) && γ >= 0) ||
        throw(ArgumentError("γ must be finite and >= 0 (photons/s), got $γ"))
    if photons !== nothing
        Base.depwarn("the photons keyword is deprecated and will be removed in 0.8.0; pass γ, the emission rate in photons/s (for new emitters γ = photons/dt gives identical output; with starting_conditions photons is ignored, as in 0.7, while γ restamps every molecule to γ·dt)", :simulate; force=true)
    end

    # Photons per record for new emitters, the rate stored in the metadata, and how it was set
    record_photons = photons !== nothing ? Float64(photons) :
                     γ !== nothing ? Float64(γ) * params.dt : 1000.0
    γ_new = photons !== nothing ? Float64(photons) / params.dt :
            γ !== nothing ? Float64(γ) : 1000.0 / params.dt
    rate_source = γ !== nothing ? "γ" : photons !== nothing ? "photons" : "default"

    # Sub-steps per exposure and per frame (warns and rounds if dt does not divide the camera timing)
    n_sub, steps_per_frame = substeps_per_frame(params)
    n_frames = floor(Int, params.t_max * params.camera_framerate + 1e-9)
    if n_frames < 1
        @warn "t_max ($(params.t_max) s) is shorter than one frame period (1/camera_framerate = $(1 / params.camera_framerate) s); simulating one frame" maxlog=1
        n_frames = 1
    end

    # Create camera if not provided
    if camera === nothing
        pixel_size = 0.1  # 100nm pixels
        n_pixels = ceil(Int, params.box_size / pixel_size)
        camera = IdealCamera(1:n_pixels, 1:n_pixels, pixel_size)
    end

    # Initialize emitters
    n_initial_emitters = 0
    prior_D = nothing
    prior_mix = nothing
    γ_restamp = γ === nothing ? nothing : Float64(γ)  # the rate every starting molecule is restamped to
    γ_val = γ_new
    if starting_conditions !== nothing
        # Extract emitters from starting_conditions
        if starting_conditions isa SMLD
            # Exact end state of the previous run, one emitter per track
            start_smld = extract_end_state(starting_conditions)
            start_emitters = start_smld.emitters
            md = start_smld.metadata
            # A saved D belongs to the tracks extract_end_state kept it for; a run without a mixture saves none
            saved_D = get(md, "monomer_D", nothing)
            prior_D = saved_D === nothing || isempty(saved_D) ? nothing : saved_D
            prior_mix = get(md, "monomer_mobility", nothing)
            if γ === nothing
                # Restamp only when the source's rate was set with γ and every molecule carries γ·dt_saved;
                # otherwise every molecule keeps its photons per record, as in 0.7.1
                saved_γ, saved_dt = get(md, "γ", nothing), get(md, "dt", nothing)
                claims_γ = get(md, "rate_source", nothing) == "γ"
                carries(e) = saved_γ isa Real && saved_dt isa Real &&
                             isapprox(e.photons, saved_γ * saved_dt; rtol=CONTINUATION_RTOL)
                if claims_γ && all(carries, start_emitters)
                    rate_source = "γ"
                    γ_restamp = Float64(saved_γ)
                else
                    claims_γ && @warn "starting_conditions: the source's rate was set with γ = $saved_γ, but not every resumed molecule carries γ·dt at its saved dt (edited, concatenated or missing dt); each keeps its photons per record, as in 0.7.1; pass γ to restamp every molecule" maxlog=1
                    rate_source = get(md, "rate_source", nothing) == "default" ? "default" : "photons"
                end
            end
        else
            # Already a vector of emitters
            start_emitters = starting_conditions
            ids = [e.track_id for e in start_emitters]
            if length(unique(ids)) != length(ids)
                @warn "starting_conditions has repeated track_id values; deduplicated to the latest record per track; pass the SMLD or extract_end_state(smld)" maxlog=1
                start_emitters = _latest_per_track(start_emitters)
            end
            γ === nothing && (rate_source = "photons")
            isempty(params.monomer_mobility) || @warn "Vector starting_conditions get fresh monomer_mobility draws; pass the SMLD from simulate (or extract_end_state(smld)) to keep each track's D" maxlog=1
        end

        # Validate emitter types
        if isempty(start_emitters)
            error("Starting conditions contain no emitters")
        end

        if !(eltype(start_emitters) <: AbstractDiffusingEmitter)
            error("Starting conditions must contain diffusing emitters")
        end

        # Every dimer needs a dimer partner that points back to it; others become monomers
        by_id = Dict(e.track_id => e for e in start_emitters)
        orphan(e) = e.state == :dimer && begin
            partner = e.partner_id === nothing ? nothing : get(by_id, e.partner_id, nothing)
            partner === nothing || partner.state != :dimer || partner.partner_id != e.track_id ||
                e.partner_id == e.track_id
        end
        if any(orphan, start_emitters)
            @warn "starting_conditions has dimers without a matching dimer partner; converted to monomers" maxlog=1
            start_emitters = [orphan(e) ? restamp(e; state=:monomer, partner_id=nothing) : e for e in start_emitters]
        end

        # Reset timestamps to start at 0.0 and frame to 1. A rate set with γ (here or by the source
        # run) restamps every molecule to γ·dt; otherwise each keeps its photons per record, as in 0.7.1
        emitters = [restamp(e; photons=(γ_restamp === nothing ? e.photons : γ_restamp * params.dt),
                            timestamp=0.0, frame=1) for e in start_emitters]
        n_initial_emitters = length(emitters)
        p1 = emitters[1].photons
        γ_val = γ_restamp !== nothing ? γ_restamp :
                all(e -> e.photons == p1, emitters) ? Float64(p1) / params.dt : nothing
    else
        # Initialize emitters using the standard approach
        emitters = build_emitters(params, record_photons, override_count)
        n_initial_emitters = length(emitters)
    end

    # Monomer diffusion coefficient per track, drawn once per molecule
    track_D = Dict{Int,Float64}()
    # A saved D is a property of the molecule and is kept even if the new mixture differs; any other track
    # draws from a non-empty mixture (saved) or uses diff_monomer at run time (never saved)
    for e in emitters
        if prior_D !== nothing && haskey(prior_D, e.track_id)
            track_D[e.track_id] = Float64(prior_D[e.track_id])
        elseif !isempty(params.monomer_mobility)
            track_D[e.track_id] = draw_monomer_D(params)
        end
    end
    if prior_D !== nothing && prior_mix !== nothing && params.monomer_mobility != prior_mix &&
       any(e -> haskey(prior_D, e.track_id), emitters)
        @warn "starting_conditions: tracks keep their saved D (metadata \"monomer_D\") although this config's monomer_mobility differs from the source run's; the new mixture applies only to tracks without a saved D; pass extract_end_state(smld).emitters for fresh draws" maxlog=1
    end

    # Store camera-frame emitters
    camera_emitters = Vector{eltype(emitters)}()

    # Simulation loop in integer steps; the first n_sub steps of each frame are recorded
    for f in 1:n_frames, j in 0:steps_per_frame-1
        k = (f - 1) * steps_per_frame + j
        j < n_sub && _record_frame!(camera_emitters, emitters, k * params.dt, f)
        emitters = update_system(emitters, params, params.dt; track_D=isempty(track_D) ? nothing : track_D)
    end

    # Convert to SMLD
    smld = create_smld(camera_emitters, camera, params; track_D=track_D, γ=γ_val, rate_source=rate_source,
                       n_frames=n_frames)

    # Live emitters are now at the start of the next frame: the exact end state for continuation
    t_end = n_frames * steps_per_frame * params.dt
    smld.metadata["last_frame_latest"] = _last_frame_latest(smld.emitters)
    smld.metadata["final_state"] = [restamp(e; timestamp=t_end, frame=n_frames) for e in emitters]

    elapsed_s = (time_ns() - start_time) / 1e9

    # Build SimInfo (diffusion doesn't have smld_true/smld_model)
    info = SimInfo(
        elapsed_s=elapsed_s,
        backend=:cpu,
        device_id=-1,
        seed=nothing,
        smld_true=nothing,
        smld_model=nothing,
        n_patterns=0,
        n_emitters=n_initial_emitters,
        n_localizations=length(smld.emitters),
        n_frames=n_frames
    )

    return smld, info
end

"""
    convert_to_diffusing_emitters(emitters::Vector{<:AbstractEmitter}, photons::Float64=1000.0, state::Symbol=:monomer)

Convert regular emitters to diffusing emitters for use as starting conditions.

# Arguments
- `emitters::Vector{<:AbstractEmitter}`: Vector of static emitters to convert
- `photons::Float64=1000.0`: Number of photons to assign (kept as given by `simulate` unless `γ` is passed)
- `state::Symbol=:monomer`: Initial state (:monomer or :dimer)

# Returns
- `Vector{<:AbstractDiffusingEmitter}`: Vector of diffusing emitters

# Example
```julia
# Convert static emitters to diffusing emitters
static_emitters = smld_static.emitters
diffusing_emitters = convert_to_diffusing_emitters(static_emitters)

# Use as starting conditions for a diffusion simulation
params = DiffusionSMLMConfig(t_max=10.0)
smld = simulate(params; starting_conditions=diffusing_emitters)
```
"""
function convert_to_diffusing_emitters(emitters::Vector{<:AbstractEmitter}, photons::Float64=1000.0, state::Symbol=:monomer)
    diffusing_emitters = Vector{Union{DiffusingEmitter2D{Float64}, DiffusingEmitter3D{Float64}}}()
    
    for (i, e) in enumerate(emitters)
        if isa(e, Emitter2D) || isa(e, Emitter2DFit)
            # Create a new diffusing emitter from the static one
            diffusing_e = DiffusingEmitter2D{Float64}(
                e.x, e.y,           # Position
                photons,            # Photons
                0.0,                # Initial timestamp
                1,                  # Initial frame
                e.dataset,          # Dataset
                e.track_id,         # ID
                state,              # Initial state
                nothing             # No partner initially
            )
            push!(diffusing_emitters, diffusing_e)
        elseif isa(e, Emitter3D) || isa(e, Emitter3DFit)
            # Create new 3D diffusing emitter
            diffusing_e = DiffusingEmitter3D{Float64}(
                e.x, e.y, e.z,      # Position
                photons,            # Photons
                0.0,                # Initial timestamp
                1,                  # Initial frame
                e.dataset,          # Dataset
                e.track_id,         # ID
                state,              # Initial state
                nothing             # No partner initially
            )
            push!(diffusing_emitters, diffusing_e)
        else
            error("Unsupported emitter type: $(typeof(e))")
        end
    end
    
    return diffusing_emitters
end

"""
    extract_end_state(smld::BasicSMLD{T,E}) where {T, E<:AbstractDiffusingEmitter}

Reduce a diffusion simulation to its exact end state, for use as `starting_conditions`.
Returns a `BasicSMLD` with one emitter per track and photons unchanged.

`simulate` stores the exact end state (the live emitters at the start of the frame after
the last) in `smld.metadata["final_state"]`, and that is returned when the SMLD's last frame is
unchanged from that run (the latest record per track of the last frame matches
`metadata["last_frame_latest"]`). Otherwise (filtered, concatenated, edited or re-wrapped
without metadata) the record with the largest timestamp per track in the last frame present is
used, with its photons as they are; tracks absent from that frame are not resumed. The result carries
`"γ"`, `"rate_source"`, `"dt"`, `"monomer_mobility"` and `"monomer_D"` when present, with `"monomer_D"`
restricted to the tracks whose resumed record is the run's own (the unchanged last frame, or a record
equal, to a relative 1e-6, to the run's stored last-frame record for that track; metadata without stored
records is taken as written), with one warning when a saved D is dropped. `simulate` then applies the
continuation rule: a γ-set rate is restamped to γ·dt only if every molecule carries γ·dt at the saved
dt, and otherwise photons per record are kept; a kept D stays, other tracks take the current setting.
Extracting twice gives the same result.

# Arguments
- `smld::BasicSMLD`: SMLD of diffusing emitters from `simulate`

# Returns
- `BasicSMLD`: One emitter per track, `n_frames = 1`, every emitter at `frame = 1`

# Example
```julia
params = DiffusionSMLMConfig(t_max=5.0)
smld, info = simulate(params)

# Continue with new parameters
params_new = DiffusionSMLMConfig(t_max=10.0, diff_monomer=0.2)
smld_continued, info = simulate(params_new; starting_conditions=extract_end_state(smld))
```
"""
function extract_end_state(smld::BasicSMLD{T,E}) where {T, E<:AbstractDiffusingEmitter}
    final = get(smld.metadata, "final_state", nothing)
    exact = final !== nothing && _final_state_matches(smld, final)
    resumed = exact ? final : _last_frame_latest(smld.emitters)
    final_emitters = [restamp(e; frame=1) for e in resumed]

    metadata = Dict{String,Any}("n_substeps" => 1)
    for key in ("simulation_type", "simulation_parameters", "dt", "camera_framerate", "camera_exposure", "γ", "rate_source",
                "monomer_mobility", "monomer_D")
        haskey(smld.metadata, key) && (metadata[key] = smld.metadata[key])
    end
    # A saved D belongs to a track only if its resumed record is the run's own
    ref = get(smld.metadata, "last_frame_latest", nothing)
    saved = get(metadata, "monomer_D", nothing)
    if !exact && ref !== nothing && saved !== nothing
        by_id = Dict(e.track_id => e for e in ref)
        own = Set(e.track_id for e in resumed if haskey(by_id, e.track_id) && _same_record(e, by_id[e.track_id]))
        dropped = count(e -> haskey(saved, e.track_id) && !(e.track_id in own), resumed)
        dropped > 0 && @warn "extract_end_state: $dropped resumed tracks with a saved D are not the run's own records (concatenated, edited or re-wrapped SMLD); they take the current diff_monomer or monomer_mobility" maxlog=1
        metadata["monomer_D"] = Dict{Int,Float64}(id => D for (id, D) in saved if id in own)
    end
    # A second extraction takes the stored path and returns the same emitters in the same order
    metadata["final_state"] = final_emitters
    metadata["last_frame_latest"] = _last_frame_latest(final_emitters)
    return BasicSMLD(final_emitters, smld.camera, 1, smld.n_datasets, metadata)
end

# Relative tolerance of the continuation checks (a molecule carries γ·dt, a record is the run's own), so an SMLD
# converted to Float32 or written and read back keeps its provenance
const CONTINUATION_RTOL = 1e-6

# Two records of one track agree: equal integer, state and partner fields, floats to CONTINUATION_RTOL. Internal.
_same_record(a, b) = nameof(typeof(a)) == nameof(typeof(b)) &&
    all(f -> _same_field(getfield(a, f), getfield(b, f)), fieldnames(typeof(a)))
_same_field(u::AbstractFloat, v::AbstractFloat) = isapprox(u, v; rtol=CONTINUATION_RTOL)
_same_field(u, v) = u == v

# The stored final state belongs to the emitters only if their last frame is the one the run
# recorded: the latest record per track there equals the stored copy. Filters, concatenation
# and edits that change the last frame fail this. Internal.
function _final_state_matches(smld::SMLD, final)
    ref = get(smld.metadata, "last_frame_latest", nothing)
    ref === nothing && return false
    return _last_frame_latest(smld.emitters) == ref
end

# One emitter per track: the record with the largest timestamp, sorted by track_id. Internal.
function _latest_per_track(emitters)
    latest = Dict{Int,eltype(emitters)}()
    for e in emitters
        if !haskey(latest, e.track_id) || e.timestamp > latest[e.track_id].timestamp
            latest[e.track_id] = e
        end
    end
    return [latest[id] for id in sort!(collect(keys(latest)))]
end

# Latest record of each track in the largest frame present. Internal.
function _last_frame_latest(emitters)
    isempty(emitters) && return similar(emitters, 0)
    max_frame = maximum(e -> e.frame, emitters)
    return _latest_per_track([e for e in emitters if e.frame == max_frame])
end

"""
    extract_final_state(smld::SMLD)

Deprecated, removed in 0.8.0: use [`extract_end_state`](@ref), which returns the exact end
state as an SMLD and carries each track's D.

Returns a `Vector`. For an SMLD of diffusing emitters, one emitter per track: the record with
the largest timestamp in the last frame, photons unchanged, sorted by track_id. For any other
SMLD (for example fitted emitters, which have no timestamp), every record in the largest frame,
as in 0.7.1.
"""
function extract_final_state(smld::SMLD)
    Base.depwarn("extract_final_state is deprecated and will be removed in 0.8.0; use extract_end_state(smld), which returns the exact end state as an SMLD and carries each track's D", :extract_final_state; force=true)
    max_frame = maximum(e -> e.frame, smld.emitters)
    return filter(e -> e.frame == max_frame, smld.emitters)
end

function extract_final_state(smld::BasicSMLD{T,E}) where {T, E<:AbstractDiffusingEmitter}
    Base.depwarn("extract_final_state is deprecated and will be removed in 0.8.0; use extract_end_state(smld), which returns the exact end state as an SMLD and carries each track's D", :extract_final_state; force=true)
    return _last_frame_latest(smld.emitters)
end
