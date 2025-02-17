# src/diffusion/smoluchowski.jl

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
end

"""
    initialize_system(params::SmoluchowskiParams)

Create initial DiffusingMoleculeSystem with randomly placed monomers.
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

"""Update position of a single monomer"""
function update_monomer_position!(mol::DiffusingMolecule, params::SmoluchowskiParams)
    σ = sqrt(2 * params.diff_monomer * params.dt)
    
    mol.x += rand(Normal(0, σ))
    mol.y += rand(Normal(0, σ))
    
    mol.updated = true
end

"""Update position and orientation of a dimer"""
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

"""Apply periodic or reflecting boundary conditions"""
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
    simulate(params::SmoluchowskiParams)

Run a complete Smoluchowski diffusion simulation.

Returns vector of DiffusingMoleculeSystem states at each timepoint.
"""
function simulate(params::SmoluchowskiParams)
    # Initialize
    n_steps = round(Int, params.t_max / params.dt)
    system = initialize_system(params)
    
    # Store system states
    systems = Vector{typeof(system)}(undef, n_steps)
    systems[1] = deepcopy(system)
    
    # Run simulation
    for t in 2:n_steps
        update_species!(system, params)
        update_positions!(system, params)
        apply_boundary!(system, params)
        
        # Store copy of current state
        systems[t] = deepcopy(system)
        systems[t].metadata["time"] = (t-1) * params.dt
    end
    
    return systems
end