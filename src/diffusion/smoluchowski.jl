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
  camera_exposure and 1/camera_framerate must be integer multiples of dt
- `t_max::Float64`: total simulation time (s); rounded down to whole frame periods
  (`t_max` shorter than one frame period throws in `simulate`)
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

Photons: the `photons` given to `simulate` is photons per emitter per frame (exposure).
Each of the `n_sub = camera_exposure/dt` records of a frame carries `photons/n_sub`.

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

# Positional construction without a mobility mixture (pre-0.8 signature)
DiffusionSMLMConfig(density, box_size, diff_monomer, diff_dimer, diff_dimer_rot,
                    k_off, r_react, d_dimer, dt, t_max, ndims, boundary,
                    camera_framerate, camera_exposure) =
    DiffusionSMLMConfig(density, box_size, diff_monomer, diff_dimer, diff_dimer_rot,
                        k_off, r_react, d_dimer, dt, t_max, ndims, boundary,
                        camera_framerate, camera_exposure, Tuple{Float64,Float64}[])

"""
    substeps_per_frame(params::DiffusionSMLMConfig) -> (n_sub, steps_per_frame)

Number of physics steps `dt` inside one exposure (`n_sub`) and in one frame period
(`steps_per_frame`). Internal. Throws `ArgumentError` unless `camera_exposure` and
`1/camera_framerate` are integer multiples of `dt` and the exposure does not exceed
the frame period.
"""
function substeps_per_frame(params::DiffusionSMLMConfig)
    r_sub = params.camera_exposure / params.dt
    r_frame = 1 / (params.camera_framerate * params.dt)
    n_sub = round(Int, r_sub)
    steps_per_frame = round(Int, r_frame)
    if abs(r_sub - n_sub) > 1e-9 * max(1, r_sub) ||
       abs(r_frame - steps_per_frame) > 1e-9 * max(1, r_frame) ||
       n_sub < 1 || n_sub > steps_per_frame
        throw(ArgumentError("dt=$(params.dt), camera_exposure=$(params.camera_exposure), " *
            "camera_framerate=$(params.camera_framerate): camera_exposure and 1/camera_framerate " *
            "must be integer multiples of dt (and exposure no longer than the frame period)"))
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
    initialize_emitters(params::DiffusionSMLMConfig, photons::Float64=1000.0; override_count::Union{Nothing, Int}=nothing)

Create initial emitter positions for the simulation.

# Arguments
- `params::DiffusionSMLMConfig`: Simulation parameters
- `photons::Float64=1000.0`: Photons per emitter per frame (exposure); each of the
  n_sub = camera_exposure/dt records in a frame carries photons/n_sub
- `override_count::Union{Nothing, Int}=nothing`: Optional override for the number of molecules

# Returns
- `Vector{<:AbstractDiffusingEmitter}`: Vector of initialized emitters
"""
function initialize_emitters(params::DiffusionSMLMConfig, photons::Float64=1000.0; override_count::Union{Nothing, Int}=nothing)
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
                photons,                   # Photons
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
                photons,                   # Photons
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
`track_D[track_id]` when given, otherwise with `params.diff_monomer`.

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
    monomer_D(id) = track_D === nothing ? params.diff_monomer : track_D[id]
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
                    error("dimer partner $(e1.partner_id) of track $(e1.track_id) not found")
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
    add_camera_frame_emitters!(camera_emitters, emitters, time, frame_num, n_sub)

Record the current emitters as one sub-step record of camera frame `frame_num`.
Each record carries `photons / n_sub`, so the `n_sub` records of a frame sum to the
per-frame photons.

# Arguments
- `camera_emitters::Vector{<:AbstractDiffusingEmitter}`: Collection of emitters for camera frames
- `emitters::Vector{<:AbstractDiffusingEmitter}`: Current emitters from simulation
- `time::Float64`: Current simulation time
- `frame_num::Int`: Current frame number
- `n_sub::Int`: Number of sub-step records per frame

# Returns
- `Nothing`
"""
function add_camera_frame_emitters!(camera_emitters, emitters, time, frame_num, n_sub)
    for e in emitters
        push!(camera_emitters, restamp(e; photons=e.photons / n_sub, timestamp=time, frame=frame_num))
    end
    
    return nothing
end

"""
    simulate(params::DiffusionSMLMConfig;
             starting_conditions::Union{Nothing, SMLD, Vector{<:AbstractDiffusingEmitter}}=nothing,
             photons::Float64=1000.0,
             override_count::Union{Nothing, Int}=nothing,
             kwargs...)

Run a Smoluchowski diffusion simulation and return a BasicSMLD object
with emitters that have both frame number and timestamp information.

# Arguments
- `params::DiffusionSMLMConfig`: Simulation parameters

# Keyword Arguments
- `starting_conditions::Union{Nothing, SMLD, Vector{<:AbstractDiffusingEmitter}}=nothing`: Optional starting emitters
  (an SMLD carries each track's D and exact end state forward; a Vector gets fresh D draws
  and must have one record per track_id)
- `photons::Float64=1000.0`: Photons per emitter per frame (exposure); each of the
  n_sub = camera_exposure/dt records in a frame carries photons/n_sub
- `override_count::Union{Nothing, Int}=nothing`: Optional override for the number of molecules
- `camera::Union{Nothing, AbstractCamera}=nothing`: Camera model (default: IdealCamera with 100nm pixels)
  - If `nothing`, creates IdealCamera with dimensions matching box_size
  - Can specify SCMOSCamera for realistic noise modeling
- Any additional parameters are ignored (allows unified interface with other simulate methods)

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
                 photons::Float64=1000.0,
                 override_count::Union{Nothing, Int}=nothing,
                 camera::Union{Nothing, AbstractCamera}=nothing,
                 kwargs...)

    start_time = time_ns()

    # Sub-steps per exposure and per frame period (validates dt against the camera timing)
    n_sub, steps_per_frame = substeps_per_frame(params)
    n_frames = floor(Int, params.t_max * params.camera_framerate + 1e-9)
    n_frames < 1 && throw(ArgumentError("t_max ($(params.t_max) s) is shorter than one frame period " *
        "(1/camera_framerate = $(1 / params.camera_framerate) s); no frames to simulate"))

    # Create camera if not provided
    if camera === nothing
        pixel_size = 0.1  # 100nm pixels
        n_pixels = ceil(Int, params.box_size / pixel_size)
        camera = IdealCamera(1:n_pixels, 1:n_pixels, pixel_size)
    end

    # Initialize emitters
    n_initial_emitters = 0
    prior_D = nothing
    if starting_conditions !== nothing
        # Extract emitters from starting_conditions
        if starting_conditions isa SMLD
            # Exact end state of the previous run, one emitter per track at per-frame photons
            start_smld = extract_final_state(starting_conditions)
            start_emitters = start_smld.emitters
            prior_D = get(start_smld.metadata, "monomer_D", nothing)
        else
            # Already a vector of emitters
            start_emitters = starting_conditions
            ids = [e.track_id for e in start_emitters]
            if length(unique(ids)) != length(ids)
                throw(ArgumentError("starting_conditions has repeated track_id values " *
                    "(one record per molecule expected); pass the SMLD from simulate or extract_final_state(smld)"))
            end
            isempty(params.monomer_mobility) || @warn "Vector starting_conditions get fresh monomer_mobility draws; pass the SMLD from simulate (or extract_final_state(smld)) to keep each track's D" maxlog=1
        end

        # Validate emitter types
        if isempty(start_emitters)
            error("Starting conditions contain no emitters")
        end

        if !(eltype(start_emitters) <: AbstractDiffusingEmitter)
            error("Starting conditions must contain diffusing emitters")
        end

        # Every dimer needs a dimer partner that points back to it
        by_id = Dict(e.track_id => e for e in start_emitters)
        for e in start_emitters
            e.state == :dimer || continue
            partner = e.partner_id === nothing ? nothing : get(by_id, e.partner_id, nothing)
            if partner === nothing || partner.state != :dimer || partner.partner_id != e.track_id
                throw(ArgumentError("starting_conditions: dimer track $(e.track_id) has no matching dimer partner " *
                    "(partner_id=$(e.partner_id))"))
            end
        end

        # Reset timestamps to start at 0.0 and frame to 1
        emitters = [restamp(e; timestamp=0.0, frame=1) for e in start_emitters]
        n_initial_emitters = length(emitters)
    else
        # Initialize emitters using the standard approach
        emitters = initialize_emitters(params, photons; override_count=override_count)
        n_initial_emitters = length(emitters)
    end

    # Monomer diffusion coefficient per track, drawn once per molecule
    track_D = Dict{Int,Float64}()
    if !isempty(params.monomer_mobility)
        for e in emitters
            track_D[e.track_id] = prior_D !== nothing && haskey(prior_D, e.track_id) ?
                Float64(prior_D[e.track_id]) : draw_monomer_D(params)
        end
    end

    # Store camera-frame emitters
    camera_emitters = Vector{eltype(emitters)}()

    # Simulation loop in integer steps; the first n_sub steps of each frame are recorded
    for f in 1:n_frames, j in 0:steps_per_frame-1
        k = (f - 1) * steps_per_frame + j
        j < n_sub && add_camera_frame_emitters!(camera_emitters, emitters, k * params.dt, f, n_sub)
        emitters = update_system(emitters, params, params.dt; track_D=isempty(track_D) ? nothing : track_D)
    end

    # Convert to SMLD
    smld = create_smld(camera_emitters, camera, params; track_D=track_D)

    # Live emitters are now at the start of the next frame: the exact end state for continuation
    t_end = n_frames * steps_per_frame * params.dt
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
- `photons::Float64=1000.0`: Number of photons to assign
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
    extract_final_state(smld::BasicSMLD{T,E}) where {T, E<:AbstractDiffusingEmitter}

Reduce a diffusion simulation to its end state, for use as `starting_conditions`.
Returns a `BasicSMLD` with one emitter per track at per-frame photons.

`simulate` stores the exact end state (the live emitters at the start of the frame after
the last) in `smld.metadata["final_state"]`, and that is returned when present. Without it
(for example an SMLD re-wrapped without metadata) the record with the largest timestamp
per track in the last frame is used, with photons multiplied by that track's number of
records in the frame. The result carries `"monomer_D"` when present, so continuation keeps
each track's D, and extracting twice gives the same result.

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
smld_continued, info = simulate(params_new; starting_conditions=extract_final_state(smld))
```
"""
function extract_final_state(smld::BasicSMLD{T,E}) where {T, E<:AbstractDiffusingEmitter}
    final = get(smld.metadata, "final_state", nothing)
    if final !== nothing
        final_emitters = [restamp(e; frame=1) for e in final]
    else
        max_frame = maximum(e -> e.frame, smld.emitters)

        # Latest record of each track in the last frame, and how many records it has there
        latest = Dict{Int,E}()
        counts = Dict{Int,Int}()
        for e in smld.emitters
            e.frame == max_frame || continue
            counts[e.track_id] = get(counts, e.track_id, 0) + 1
            if !haskey(latest, e.track_id) || e.timestamp > latest[e.track_id].timestamp
                latest[e.track_id] = e
            end
        end

        final_emitters = [restamp(latest[id]; photons=latest[id].photons * counts[id], frame=1)
                          for id in sort!(collect(keys(latest)))]
    end

    metadata = Dict{String,Any}("n_substeps" => 1, "final_state" => final_emitters)
    for key in ("simulation_type", "simulation_parameters", "camera_framerate", "camera_exposure", "monomer_D")
        haskey(smld.metadata, key) && (metadata[key] = smld.metadata[key])
    end
    return BasicSMLD(final_emitters, smld.camera, 1, smld.n_datasets, metadata)
end
