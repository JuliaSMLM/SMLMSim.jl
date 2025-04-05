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
Base.@kwdef mutable struct SmoluchowskiParams
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
    initialize_system(params::SmoluchowskiParams)

Create initial DiffusingMoleculeSystem with randomly placed monomers.

# Arguments
- `params::SmoluchowskiParams`: Simulation parameters

# Returns
- `DiffusingMoleculeSystem`: Initialized system with molecules
"""
function initialize_system(params::SmoluchowskiParams)
    # Calculate number of molecules
    n_molecules = round(Int, params.density * params.box_size^params.ndims)
    
    # Create camera matching box size
    pixel_size = 0.1  # 100nm pixels
    n_pixels = ceil(Int, params.box_size / pixel_size)
    camera = IdealCamera(1:n_pixels, 1:n_pixels, pixel_size)
    
    # Create emitters at random positions
    molecules = Vector{DiffusingMolecule{Emitter2D{Float64}}}(undef, n_molecules)
    for i in 1:n_molecules
        x = rand(Uniform(0, params.box_size))
        y = rand(Uniform(0, params.box_size))
        
        # Create basic emitter with position and standard brightness
        emitter = Emitter2D{Float64}(x, y, 1000.0)
        molecules[i] = DiffusingMolecule(emitter, 1, i, nothing, false)
    end
    
    # Create system
    DiffusingMoleculeSystem(
        molecules,
        camera,
        params.box_size,
        1,  # Start with single frame
        1,  # Single dataset
        Dict{String,Any}(
            "simulation_parameters" => params
        )
    )
end

"""
    update_species!(system::DiffusingMoleculeSystem, params::SmoluchowskiParams)

Update molecular states (dimerization/dissociation) for all molecules.

# Arguments
- `system::DiffusingMoleculeSystem`: System to update
- `params::SmoluchowskiParams`: Simulation parameters

# Returns
- `Nothing`

# Details
For each molecule:
1. If monomer, check for possible dimerization with other monomers
2. If dimer, check for possible dissociation based on k_off rate
"""
function update_species!(system::DiffusingMoleculeSystem, params::SmoluchowskiParams)
    for (i, mol1) in enumerate(system.molecules)
        if mol1.state == 1  # monomer
            # Check for dimerization with other monomers
            for mol2 in @view system.molecules[i+1:end]
                if mol2.state == 1  # also a monomer
                    if calc_r(mol1, mol2) < params.r_react
                        dimerize!(mol1, mol2, params.d_dimer)
                        break  # Only one dimerization per molecule per step
                    end
                end
            end
        elseif mol1.state == 2  # dimer
            # Check for dissociation
            if rand() < params.k_off * params.dt
                monomerize!(mol1, system)
            end
        end
    end
end

"""
    update_positions!(system::DiffusingMoleculeSystem, params::SmoluchowskiParams)

Update positions of all molecules using appropriate diffusion models.

# Arguments
- `system::DiffusingMoleculeSystem`: System to update
- `params::SmoluchowskiParams`: Simulation parameters

# Returns
- `Nothing`

# Details
For each molecule:
1. If monomer, update position with Brownian motion
2. If dimer, update position and orientation with coupled diffusion
3. Molecule positions are only updated once per timestep
"""
function update_positions!(system::DiffusingMoleculeSystem, params::SmoluchowskiParams)
    for mol in system.molecules
        if !mol.updated  # Skip if already updated through dimer
            if mol.state == 1
                update_monomer_position!(mol, params)
            else  # state == 2
                # Find linked molecule and update dimer
                linked_idx = findfirst(m -> m.id == mol.link, system.molecules)
                if !isnothing(linked_idx)
                    update_dimer_position!(mol, system.molecules[linked_idx], params)
                end
            end
        end
    end
    
    # Reset update flags
    for mol in system.molecules
        mol.updated = false
    end
end

"""
    update_monomer_position!(mol::DiffusingMolecule, params::SmoluchowskiParams)

Update position of a single monomer using Brownian motion.

# Arguments
- `mol::DiffusingMolecule`: Molecule to update
- `params::SmoluchowskiParams`: Simulation parameters

# Returns
- `Nothing`
"""
function update_monomer_position!(mol::DiffusingMolecule, params::SmoluchowskiParams)
    σ = sqrt(2 * params.diff_monomer * params.dt)
    
    mol.x += rand(Normal(0, σ))
    mol.y += rand(Normal(0, σ))
    
    mol.updated = true
end

"""
    update_dimer_position!(mol1::DiffusingMolecule, mol2::DiffusingMolecule, params::SmoluchowskiParams)

Update position and orientation of a dimer with both translational and rotational diffusion.

# Arguments
- `mol1::DiffusingMolecule`: First molecule in dimer
- `mol2::DiffusingMolecule`: Second molecule in dimer
- `params::SmoluchowskiParams`: Simulation parameters

# Returns
- `Nothing`

# Details
1. Translational diffusion of center of mass
2. Rotational diffusion around center of mass
3. Maintains fixed separation distance between molecules
"""
function update_dimer_position!(mol1::DiffusingMolecule, mol2::DiffusingMolecule, params::SmoluchowskiParams)
    # Translational diffusion of center of mass
    σ_trans = sqrt(2 * params.diff_dimer * params.dt)
    dx = rand(Normal(0, σ_trans))
    dy = rand(Normal(0, σ_trans))
    
    # Current center of mass
    com_x = (mol1.x + mol2.x)/2
    com_y = (mol1.y + mol2.y)/2
    
    # Move center of mass
    com_x += dx
    com_y += dy
    
    # Rotational diffusion
    σ_rot = sqrt(2 * params.diff_dimer_rot * params.dt)
    ϕ = calc_ϕ(mol1, mol2)
    ϕ += rand(Normal(0, σ_rot))
    
    # Calculate new positions relative to center of mass
    r = params.d_dimer/2
    dx = r * cos(ϕ)
    dy = r * sin(ϕ)
    
    # Update positions
    mol1.x = com_x - dx
    mol1.y = com_y - dy
    mol2.x = com_x + dx
    mol2.y = com_y + dy
    
    mol1.updated = true
    mol2.updated = true
end

"""
    apply_boundary!(system::DiffusingMoleculeSystem, params::SmoluchowskiParams)

Apply periodic or reflecting boundary conditions to keep molecules in simulation box.

# Arguments
- `system::DiffusingMoleculeSystem`: System to apply boundaries to
- `params::SmoluchowskiParams`: Simulation parameters

# Returns
- `Nothing`

# Details
- For periodic boundaries, positions wrap around the box edges
- For reflecting boundaries, positions are reflected back into the box
"""
function apply_boundary!(system::DiffusingMoleculeSystem, params::SmoluchowskiParams)
    for mol in system.molecules
        if params.boundary == "periodic"
            mol.x = mod(mol.x, params.box_size)
            mol.y = mod(mol.y, params.box_size)
        else  # reflecting
            if mol.x < 0
                mol.x = -mol.x
            elseif mol.x > params.box_size
                mol.x = 2params.box_size - mol.x
            end
            
            if mol.y < 0
                mol.y = -mol.y
            elseif mol.y > params.box_size
                mol.y = 2params.box_size - mol.y
            end
        end
    end
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
    system = initialize_system(params)
    
    # Use appropriate emitter type based on dimensions
    EmitterType = params.ndims == 3 ? DiffusingEmitter3D{Float64} : DiffusingEmitter2D{Float64}
    all_emitters = Vector{EmitterType}()
    
    # Calculate maximum frame number
    max_frame = ceil(Int, params.t_max * params.camera_framerate)
    
    # Add initial state to emitters - need to modify for 3D if ndims=3
    time_val = 0.0
    frame_num = 1
    
    # Check if initial timepoint falls within a camera exposure window
    exposure_start = 0.0
    exposure_end = params.camera_exposure
    
    if time_val >= exposure_start && time_val <= exposure_end
        for mol in system.molecules
            if params.ndims == 2
                # Create 2D emitter
                emitter = DiffusingEmitter2D{Float64}(
                    mol.x, mol.y,        # spatial coordinates
                    photons,             # photon count
                    time_val,            # actual timestamp
                    frame_num,           # frame number
                    1,                   # dataset
                    mol.id,              # track_id 
                    mol.state,           # monomer/dimer state
                    mol.link             # linked molecule id
                )
            else
                # Create 3D emitter (assuming z=0 for now, can be extended for 3D)
                emitter = DiffusingEmitter3D{Float64}(
                    mol.x, mol.y, 0.0,   # spatial coordinates
                    photons,             # photon count
                    time_val,            # actual timestamp
                    frame_num,           # frame number
                    1,                   # dataset
                    mol.id,              # track_id 
                    mol.state,           # monomer/dimer state
                    mol.link             # linked molecule id
                )
            end
            push!(all_emitters, emitter)
        end
    end
    
    # Run simulation
    for t in 2:n_steps
        # Update molecule positions and states
        update_species!(system, params)
        update_positions!(system, params)
        apply_boundary!(system, params)
        
        # Calculate time and frame for this timepoint
        time_val = (t-1) * params.dt
        frame_num = ceil(Int, time_val * params.camera_framerate)
        
        # Skip if frame is past max frames
        frame_num > max_frame && continue
        
        # Check if this timepoint falls within a camera exposure window
        exposure_start = (frame_num - 1) / params.camera_framerate
        exposure_end = exposure_start + params.camera_exposure
        
        if time_val >= exposure_start && time_val <= exposure_end
            # Create emitters for all molecules at this timepoint
            for mol in system.molecules
                if params.ndims == 2
                    # Create 2D emitter
                    emitter = DiffusingEmitter2D{Float64}(
                        mol.x, mol.y,        # spatial coordinates
                        photons,             # photon count
                        time_val,            # actual timestamp
                        frame_num,           # frame number
                        1,                   # dataset
                        mol.id,              # track_id 
                        mol.state,           # monomer/dimer state
                        mol.link             # linked molecule id
                    )
                else
                    # Create 3D emitter (assuming z=0 for now, can be extended for 3D)
                    emitter = DiffusingEmitter3D{Float64}(
                        mol.x, mol.y, 0.0,   # spatial coordinates
                        photons,             # photon count
                        time_val,            # actual timestamp
                        frame_num,           # frame number
                        1,                   # dataset
                        mol.id,              # track_id 
                        mol.state,           # monomer/dimer state
                        mol.link             # linked molecule id
                    )
                end
                push!(all_emitters, emitter)
            end
        end
    end
    
    # Create a single SMLD with all emitters
    return BasicSMLD(
        all_emitters,
        system.camera,
        max_frame,       # nframes based on camera framerate
        1,               # ndatasets
        Dict{String,Any}(
            "simulation_parameters" => params,
            "simulation_type" => "diffusion",
            "camera_framerate" => params.camera_framerate,
            "camera_exposure" => params.camera_exposure
        )
    )
end