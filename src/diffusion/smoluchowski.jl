"""
    SmoluchowskiParams

Parameters for Smoluchowski diffusion simulation.

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
params = SmoluchowskiParams()

# Custom parameters
params = SmoluchowskiParams(
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
Base.@kwdef mutable struct SmoluchowskiParams <: AbstractSim
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
    
    function SmoluchowskiParams(
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
    initialize_emitters(params::SmoluchowskiParams, photons::Float64=1000.0)

Create initial emitter positions for the simulation.

# Arguments
- `params::SmoluchowskiParams`: Simulation parameters
- `photons::Float64=1000.0`: Number of photons per emitter

# Returns
- `Vector{<:AbstractDiffusingEmitter}`: Vector of initialized emitters
"""
function initialize_emitters(params::SmoluchowskiParams, photons::Float64=1000.0)
    # Calculate number of molecules
    n_molecules = round(Int, params.density * params.box_size^params.ndims)
    
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
                i,                         # ID
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
                i,                         # ID
                :monomer,                  # Initial state
                nothing                    # No partner initially
            )
        end
    end
    
    return emitters
end

"""
    update_system(emitters::Vector{<:AbstractDiffusingEmitter}, params::SmoluchowskiParams, dt::Float64) 

Update all emitters based on Smoluchowski diffusion dynamics.

# Arguments
- `emitters::Vector{<:AbstractDiffusingEmitter}`: Current emitters state
- `params::SmoluchowskiParams`: Simulation parameters
- `dt::Float64`: Time step

# Returns
- `Vector{<:AbstractDiffusingEmitter}`: Updated emitters
"""
function update_system(emitters::Vector{<:AbstractDiffusingEmitter}, params::SmoluchowskiParams, dt::Float64)
    # Create new array for updated emitters
    new_emitters = Vector{eltype(emitters)}()
    
    # Create a set of IDs that have been processed
    processed = Set{Int}()
    
    # Process all emitters
    for (i, e1) in enumerate(emitters)
        e1.id in processed && continue
        
        if e1.state == :monomer
            # Check if this monomer forms a dimer with any other monomer
            found_dimer = false
            
            for (j, e2) in enumerate(emitters[i+1:end])
                e2.id in processed && continue
                e2.state == :monomer || continue
                
                if can_dimerize(e1, e2, params.r_react)
                    # Create new dimer pair
                    d1, d2 = dimerize(e1, e2, params.d_dimer)
                    push!(new_emitters, d1, d2)
                    push!(processed, e1.id, e2.id)
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
                push!(processed, e1.id)
            end
        elseif e1.state == :dimer && !(e1.id in processed)
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
                push!(processed, e1.id, e1.partner_id)
            else
                # Find partner and update dimer
                partner_idx = findfirst(e -> e.id == e1.partner_id, emitters)
                if !isnothing(partner_idx) && !(emitters[partner_idx].id in processed)
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
                    push!(processed, e1.id, e2.id)
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
- `params::SmoluchowskiParams`: Simulation parameters

# Returns
- `Nothing`
"""
function add_camera_frame_emitters!(camera_emitters, emitters, time, frame_num, params)
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
                    e.id,               # ID
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
                    e.id,               # ID
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
    simulate(params::SmoluchowskiParams; photons::Float64=1000.0, kwargs...)

Run a Smoluchowski diffusion simulation and return a BasicSMLD object
with emitters that have both frame number and timestamp information.

# Arguments
- `params::SmoluchowskiParams`: Simulation parameters
- `photons::Float64=1000.0`: Number of photons per emitter

# Keyword Arguments
- Any additional parameters are ignored (allows unified interface with other simulate methods)

# Returns
- `BasicSMLD`: Single SMLD object containing all emitters across all frames

# Example
```julia
# Set up parameters with camera settings
params = SmoluchowskiParams(
    density = 0.5,           # molecules per μm²
    box_size = 10.0,         # μm
    camera_framerate = 20.0, # 20 fps
    camera_exposure = 0.04   # 40ms exposure
)

# Run simulation
smld = simulate(params)

# Generate images
psf = Gaussian2D(0.15)  # 150nm PSF width
images = gen_images(psf, smld)
```
"""
function simulate(params::SmoluchowskiParams; photons::Float64=1000.0, kwargs...)
    # Initialize
    n_steps = round(Int, params.t_max / params.dt)
    
    # Create camera
    pixel_size = 0.1  # 100nm pixels
    n_pixels = ceil(Int, params.box_size / pixel_size)
    camera = IdealCamera(1:n_pixels, 1:n_pixels, pixel_size)
    
    # Initialize emitters
    emitters = initialize_emitters(params, photons)
    
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
