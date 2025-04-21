"""
    DiffusionSMLMParams <: SMLMSimParams

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
- `dt::Float64`: time step (s)
- `t_max::Float64`: total simulation time (s)
- `ndims::Int`: number of dimensions (2 or 3)
- `boundary::String`: boundary condition type ("periodic" or "reflecting")
- `camera_framerate::Float64`: camera frames per second (Hz)
- `camera_exposure::Float64`: camera exposure time per frame (s)

# Examples
```julia
# Default parameters
params = DiffusionSMLMParams()

# Custom parameters
params = DiffusionSMLMParams(
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
```
"""
Base.@kwdef mutable struct DiffusionSMLMParams <: SMLMSimParams
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
    
    function DiffusionSMLMParams(
        density, box_size, diff_monomer, diff_dimer, diff_dimer_rot,
        k_off, r_react, d_dimer, dt, t_max, ndims, boundary,
        camera_framerate, camera_exposure
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
        
        new(density, box_size, diff_monomer, diff_dimer, diff_dimer_rot,
            k_off, r_react, d_dimer, dt, t_max, ndims, boundary,
            camera_framerate, camera_exposure)
    end
end

"""
    initialize_emitters(params::DiffusionSMLMParams, photons::Float64=1000.0; override_count::Union{Nothing, Int}=nothing)

Create initial emitter positions for the simulation.

# Arguments
- `params::DiffusionSMLMParams`: Simulation parameters
- `photons::Float64=1000.0`: Number of photons per emitter
- `override_count::Union{Nothing, Int}=nothing`: Optional override for the number of molecules

# Returns
- `Vector{<:AbstractDiffusingEmitter}`: Vector of initialized emitters
"""
function initialize_emitters(params::DiffusionSMLMParams, photons::Float64=1000.0; override_count::Union{Nothing, Int}=nothing)
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
    update_system(emitters::Vector{<:AbstractDiffusingEmitter}, params::DiffusionSMLMParams, dt::Float64) 

Update all emitters based on Smoluchowski diffusion dynamics.

# Arguments
- `emitters::Vector{<:AbstractDiffusingEmitter}`: Current emitters state
- `params::DiffusionSMLMParams`: Simulation parameters
- `dt::Float64`: Time step

# Returns
- `Vector{<:AbstractDiffusingEmitter}`: Updated emitters
"""
function update_system(emitters::Vector{<:AbstractDiffusingEmitter}, params::DiffusionSMLMParams, dt::Float64)
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
                new_e = diffuse(e1, params.diff_monomer, dt)
                
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
                m1 = diffuse(m1, params.diff_monomer, dt)
                m2 = diffuse(m2, params.diff_monomer, dt)
                
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
                end
            end
        end
    end
    
    return new_emitters
end

"""
    add_camera_frame_emitters!(camera_emitters, emitters, time, frame_num, params)

Add emitters to camera frames when they fall within an exposure window.

# Arguments
- `camera_emitters::Vector{<:AbstractDiffusingEmitter}`: Collection of emitters for camera frames
- `emitters::Vector{<:AbstractDiffusingEmitter}`: Current emitters from simulation
- `time::Float64`: Current simulation time
- `frame_num::Int`: Current frame number
- `params::DiffusionSMLMParams`: Simulation parameters

# Returns
- `Nothing`
"""
function add_camera_frame_emitters!(camera_emitters, emitters, time, frame_num, params::DiffusionSMLMParams)
    # Check if this timepoint falls within a camera exposure window
    exposure_start = (frame_num - 1) / params.camera_framerate
    exposure_end = exposure_start + params.camera_exposure
    
    if time >= exposure_start && time <= exposure_end
        # Add all current emitters to the camera frame
        for e in emitters
            # Create a copy with the correct frame number
            if isa(e, DiffusingEmitter2D)
                camera_emitter = DiffusingEmitter2D{typeof(e.x)}(
                    e.x, e.y,           # Position
                    e.photons,          # Photons
                    time,               # Current timestamp
                    frame_num,          # Frame number
                    e.dataset,          # Dataset
                    e.track_id,               # ID
                    e.state,            # State
                    e.partner_id        # Partner ID
                )
            else  # 3D
                camera_emitter = DiffusingEmitter3D{typeof(e.x)}(
                    e.x, e.y, e.z,      # Position
                    e.photons,          # Photons
                    time,               # Current timestamp
                    frame_num,          # Frame number
                    e.dataset,          # Dataset
                    e.track_id,               # ID
                    e.state,            # State
                    e.partner_id        # Partner ID
                )
            end
            push!(camera_emitters, camera_emitter)
        end
    end
    
    return nothing
end

"""
    simulate(params::DiffusionSMLMParams; 
             starting_conditions::Union{Nothing, SMLD, Vector{<:AbstractDiffusingEmitter}}=nothing,
             photons::Float64=1000.0, 
             override_count::Union{Nothing, Int}=nothing,
             kwargs...)

Run a Smoluchowski diffusion simulation and return a BasicSMLD object
with emitters that have both frame number and timestamp information.

# Arguments
- `params::DiffusionSMLMParams`: Simulation parameters
- `starting_conditions::Union{Nothing, SMLD, Vector{<:AbstractDiffusingEmitter}}`: Optional starting emitters
- `photons::Float64=1000.0`: Number of photons per emitter
- `override_count::Union{Nothing, Int}=nothing`: Optional override for the number of molecules

# Keyword Arguments
- Any additional parameters are ignored (allows unified interface with other simulate methods)

# Returns
- `BasicSMLD`: Single SMLD object containing all emitters across all frames

# Example
```julia
# Set up parameters with camera settings
params = DiffusionSMLMParams(
    density = 0.5,           # molecules per μm²
    box_size = 10.0,         # μm
    camera_framerate = 20.0, # 20 fps
    camera_exposure = 0.04   # 40ms exposure
)

# Run basic simulation
smld = simulate(params)

# Run simulation with exactly 2 particles
smld = simulate(params; override_count=2)

# Use previous simulation state as starting conditions for a new simulation
final_frame = maximum([e.frame for e in smld.emitters])
final_state_emitters = filter(e -> e.frame == final_frame, smld.emitters)
smld_continued = simulate(params; starting_conditions=final_state_emitters)
```
"""
function simulate(params::DiffusionSMLMParams; 
                 starting_conditions::Union{Nothing, SMLD, Vector{<:AbstractDiffusingEmitter}}=nothing,
                 photons::Float64=1000.0, 
                 override_count::Union{Nothing, Int}=nothing,
                 kwargs...)
    # Initialize
    n_steps = round(Int, params.t_max / params.dt)
    
    # Create camera
    pixel_size = 0.1  # 100nm pixels
    n_pixels = ceil(Int, params.box_size / pixel_size)
    camera = IdealCamera(1:n_pixels, 1:n_pixels, pixel_size)
    
    # Initialize emitters
    if starting_conditions !== nothing
        # Extract emitters from starting_conditions
        if starting_conditions isa SMLD
            # Get emitters from SMLD, using the most recent frame
            max_frame = maximum([e.frame for e in starting_conditions.emitters])
            start_emitters = filter(e -> e.frame == max_frame, starting_conditions.emitters)
        else
            # Already a vector of emitters
            start_emitters = starting_conditions
        end
        
        # Validate emitter types
        if isempty(start_emitters)
            error("Starting conditions contain no emitters")
        end
        
        if !(eltype(start_emitters) <: AbstractDiffusingEmitter)
            error("Starting conditions must contain diffusing emitters")
        end
        
        # Create deep copies of the starting emitters
        emitters = deepcopy.(start_emitters)
        
        # Reset timestamps to start at 0.0 and frame to 1
        for i in eachindex(emitters)
            if isa(emitters[i], DiffusingEmitter2D)
                emitters[i] = DiffusingEmitter2D{typeof(emitters[i].x)}(
                    emitters[i].x, emitters[i].y,  # Position
                    emitters[i].photons,           # Photons
                    0.0,                           # Reset timestamp to 0
                    1,                             # Initial frame
                    emitters[i].dataset,           # Dataset
                    emitters[i].track_id,                # ID
                    emitters[i].state,             # State
                    emitters[i].partner_id         # Partner ID
                )
            elseif isa(emitters[i], DiffusingEmitter3D)
                emitters[i] = DiffusingEmitter3D{typeof(emitters[i].x)}(
                    emitters[i].x, emitters[i].y, emitters[i].z,  # Position
                    emitters[i].photons,                         # Photons
                    0.0,                                         # Reset timestamp to 0
                    1,                                           # Initial frame
                    emitters[i].dataset,                         # Dataset
                    emitters[i].track_id,                              # ID
                    emitters[i].state,                           # State
                    emitters[i].partner_id                       # Partner ID
                )
            end
        end
    else
        # Initialize emitters using the standard approach
        emitters = initialize_emitters(params, photons; override_count=override_count)
    end
    
    # Store camera-frame emitters
    camera_emitters = Vector{eltype(emitters)}()
    
    # Add initial state to camera frames if in exposure window
    add_camera_frame_emitters!(camera_emitters, emitters, 0.0, 1, params)
    
    # Simulation loop
    time = 0.0
    while time < params.t_max
        # Update time
        time += params.dt
        
        # Update emitter states and positions
        emitters = update_system(emitters, params, params.dt)
        
        # Calculate frame number for this timepoint
        frame_num = ceil(Int, time * params.camera_framerate)
        
        # Add to camera frames if in exposure window
        add_camera_frame_emitters!(camera_emitters, emitters, time, frame_num, params)
    end
    
    # Convert to SMLD
    return create_smld(camera_emitters, camera, params)
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
params = DiffusionSMLMParams(t_max=10.0)
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
    extract_final_state(smld::SMLD)

Extract the emitters from the final frame of a simulation to use as starting conditions.

# Arguments
- `smld::SMLD`: SMLD containing emitters from a simulation

# Returns
- `Vector{<:AbstractEmitter}`: Emitters from the final frame

# Example
```julia
# Run a simulation
params = DiffusionSMLMParams(t_max=5.0)
smld = simulate(params)

# Extract final state
final_state = extract_final_state(smld)

# Continue simulation with new parameters
params_new = DiffusionSMLMParams(t_max=10.0, diff_monomer=0.2)
smld_continued = simulate(params_new; starting_conditions=final_state)
```
"""
function extract_final_state(smld::SMLD)
    # Find the maximum frame number
    max_frame = maximum([e.frame for e in smld.emitters])
    
    # Filter emitters from the maximum frame
    final_emitters = filter(e -> e.frame == max_frame, smld.emitters)
    
    return final_emitters
end
