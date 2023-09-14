
abstract type AbstractOligomer end

# kwargs struct for Smoluchowski simulation
Base.@kwdef mutable struct ArgsSmol
    density::Float64 = 1.0
    box_size::Float64 = 10.0
    diff_monomer::Float64 = 0.1
    diff_dimer::Float64 = 0.05
    diff_dimer_rot::Union{Nothing,Float64} = nothing
    k_off::Float64 = 0.2
    r_react::Float64 = 0.01
    dt::Float64 = 0.01
    t_max::Float64 = 100.0
    ndims::Int64 = 2
    d_dimer::Float64 = 0.01
    boundary::String = "periodic"
end


"""
    struct Monomer <: AbstractOligomer

A struct representing a single monomer in an oligomer.

# Fields
- `x::Float64`: the x-coordinate of the monomer's position
- `y::Float64`: the y-coordinate of the monomer's position
- `z::Float64`: the z-coordinate of the monomer's position
- `state::Int64`: the state of the monomer (1 for monomer, 2 for dimer)
- `link::Union{Monomer,Nothing}`: the monomer that this monomer is linked to, if any
- `updated::Bool`: a flag indicating whether this monomer has been updated in the current iteration
"""
mutable struct Monomer <: AbstractOligomer
    x::Float64
    y::Float64
    z::Float64
    state::Int64
    link::Union{Monomer,Nothing}
    updated::Bool
end

"""
    struct MoleculeStates

A struct representing the states of a system of molecules.

# Fields
- `dt::Float64`: the time step used in the simulation
- `States::Vector{Vector{Monomer}}`: a vector of vectors of `Monomer` objects representing the states of the molecules in the system at each time step
"""
struct MoleculeStates
    dt::Float64
    States::Vector{Vector{Monomer}}
end

"""
    dimerize!(mol1::Monomer, mol2::Monomer, distance::Float64)

Update the state and position of two monomers to form a dimer.

# Arguments
- `mol1::Monomer`: the first monomer to dimerize
- `mol2::Monomer`: the second monomer to dimerize
- `distance::Float64`: the distance between the monomers in the dimer

# Returns
- `nothing`

"""
function dimerize!(mol1::Monomer, mol2::Monomer, distance::Float64)

    # Update status
    mol1.state = 2
    mol2.state = 2
    mol1.link = mol2
    mol2.link = mol1

    # Calculate center of mass
    com_x = (mol1.x + mol2.x) / 2
    com_y = (mol1.y + mol2.y) / 2
    com_z = (mol1.z + mol2.z) / 2

    # Find angle between monomers
    ϕ = calc_ϕ(mol1, mol2)
    θ = calc_θ(mol1, mol2)
    r = distance / 2

    # Set positions of monomers using center of mass and angle and distance
    x = r * cos(ϕ) * sin(θ)
    y = r * sin(ϕ) * sin(θ)
    z = r * cos(θ)

    # Update positions
    mol1.x = com_x - x
    mol1.y = com_y - y
    mol1.z = com_z - z
    mol2.x = com_x + x
    mol2.y = com_y + y
    mol2.z = com_z + z

    return nothing
end

"""
    monomerize!(mol::Monomer)

Update the state and links of a monomer and its linked molecule.

# Arguments
- `mol::Monomer`: the monomer to monomerize

# Returns
- `nothing`

"""
function monomerize!(mol::Monomer)
    mol.state = 1
    mol.link.state = 1
    mol.link.link = nothing
    mol.link = nothing
    return nothing
end


"""
    calc_r(mol1::Monomer, mol2::Monomer)

Calculate the Euclidean distance between two monomers.

# Arguments
- `mol1::Monomer`: the first monomer
- `mol2::Monomer`: the second monomer

# Returns
- `r::Float64`: the Euclidean distance between the two monomers

"""
function calc_r(mol1::Monomer, mol2::Monomer)
    return sqrt((mol1.x - mol2.x)^2 + (mol1.y - mol2.y)^2 + (mol1.z - mol2.z)^2)
end

"""
    calc_θ(mol1::Monomer, mol2::Monomer)

Calculate the angle between the z-axis and the vector connecting two monomers.

# Arguments
- `mol1::Monomer`: the first monomer
- `mol2::Monomer`: the second monomer

# Returns
- `θ::Float64`: the angle between the z-axis and the vector connecting the two monomers, in radians

"""
function calc_θ(mol1::Monomer, mol2::Monomer)
    x = mol2.x - mol1.x
    y = mol2.y - mol1.y
    z = mol2.z - mol1.z
    r = sqrt(x^2 + y^2 + z^2)
    return acos(z / r)
end

"""
    calc_ϕ(mol1::Monomer, mol2::Monomer)

Calculate the azimuthal angle between two monomers.

# Arguments
- `mol1::Monomer`: the first monomer
- `mol2::Monomer`: the second monomer

# Returns
- `ϕ::Float64`: the azimuthal angle between the two monomers, in radians

"""
function calc_ϕ(mol1::Monomer, mol2::Monomer)
    # Calculate azimuthal angle between monomers
    x = mol2.x - mol1.x
    y = mol2.y - mol1.y
    return atan(y, x)
end

"""
    update_species!(molecules::Vector{<:AbstractOligomer}, args::ArgsSmol)

Update the state of a system of molecules according to the Smoluchowski model.

# Arguments
- `molecules::Vector{<:AbstractOligomer}`: a vector of `AbstractOligomer` objects representing the molecules in the system
- `args::ArgsSmol`: a struct containing the simulation parameters

# Returns
- `nothing`

"""
function update_species!(molecules::Vector{<:AbstractOligomer}, args::ArgsSmol)
    for mol1 in molecules
        if mol1.state == 1 # a monomer 
            # Check to see if any other molecules are within r_react
            # Only check for molecules that are later in the list
            # to avoid double counting
            for mol2 in molecules[findfirst(x -> x == mol1, molecules)+1:end]
                if mol2.state == 1 # a monomer
                    if (mol1.x - mol2.x)^2 + (mol1.y - mol2.y)^2 + (mol1.z - mol2.z)^2 < args.r_react^2
                        # Form dimer by updating links and enforcing distance
                        dimerize!(mol1, mol2, args.d_dimer)
                    end
                end
            end

        elseif mol1.state == 2 # a dimer
            if rand() < args.k_off * args.dt
                monomerize!(mol1) # converts both monomers to state 1
            end
        end
    end
    return nothing
end

"""
    update_position!(mol::Monomer, args::ArgsSmol)

Update the position of a single monomer using isotropic Brownian motion.

# Arguments
- `mol::Monomer`: the monomer to update
- `args::ArgsSmol`: a struct containing the simulation parameters

# Returns
- `nothing`

"""
function update_position!(mol::Monomer, args::ArgsSmol)
    # Update position of a single molecule
    diff = args.diff_monomer
    dt = args.dt
    mol.x += rand(Normal(0.0, sqrt(2 * diff * dt)))
    mol.y += rand(Normal(0.0, sqrt(2 * diff * dt)))
    if args.ndims == 3
        mol.z += rand(Normal(0.0, sqrt(2 * diff * dt)))
    end
    mol.updated = true
    return nothing
end

"""
    update_position!(mol1::Monomer, mol2::Monomer, args::ArgsSmol)

Update the position of a dimer using isotropic Brownian motion.

The center of mass of the dimer is updated using isotropic Brownian motion,
and the orientation of the dimer is updated using rotational Brownian motion.

# Arguments
- `mol1::Monomer`: the first monomer in the dimer to update
- `mol2::Monomer`: the second monomer in the dimer to update
- `args::ArgsSmol`: a struct containing the simulation parameters

# Returns
- `nothing`

"""
function update_position!(mol1::Monomer, mol2::Monomer, args::ArgsSmol)
    # Update center of mass
    diff = args.diff_dimer
    dt = args.dt
    com_x = (mol1.x + mol2.x) / 2
    com_y = (mol1.y + mol2.y) / 2
    com_z = (mol1.z + mol2.z) / 2
    com_x += rand(Normal(0.0, sqrt(2 * diff * dt)))
    com_y += rand(Normal(0.0, sqrt(2 * diff * dt)))
    if args.ndims == 3
        com_z += rand(Normal(0.0, sqrt(2 * diff * dt)))
    end

    # Update orientation
    diff_rot = args.diff_dimer_rot
    
    ϕ = calc_ϕ(mol1, mol2)
    ϕ += rand(Normal(0.0, sqrt(2 * diff_rot * dt)))

    if args.ndims == 3 # 3D
        θ = calc_θ(mol1, mol2)
        θ += rand(Normal(0.0, sqrt(2 * diff_rot * dt)))
        mol1.z = com_z - args.d_dimer * sin(θ) * cos(ϕ)
        mol2.z = com_z + args.d_dimer * sin(θ) * cos(ϕ)
    end

    mol1.x = com_x - args.d_dimer * cos(ϕ)
    mol1.y = com_y - args.d_dimer * sin(ϕ)
    mol2.x = com_x + args.d_dimer * cos(ϕ)
    mol2.y = com_y + args.d_dimer * sin(ϕ)

    mol1.updated = true
    mol2.updated = true
    return nothing
end

"""
    update_positions!(molecules::Vector{<:AbstractOligomer}, args::ArgsSmol)

Update the positions of all molecules in the system according to the Smoluchowski model.

# Arguments
- `molecules::Vector{<:AbstractOligomer}`: a vector of `AbstractOligomer` objects representing the molecules in the system
- `args::ArgsSmol`: a struct containing the simulation parameters

# Returns
- `nothing`

"""
function update_positions!(molecules::Vector{<:AbstractOligomer}, args::ArgsSmol)
    # Update positions of all molecules
    for mol in molecules
        if mol.state == 1 && mol.updated == false
            update_position!(mol, args)
        elseif mol.state == 2
            update_position!(mol, mol.link, args)
        end
    end

    for mol in molecules
        mol.updated = false
    end

    return nothing
end

"""
    apply_boundary!(molecules::Vector{<:AbstractOligomer}, args::ArgsSmol)

Apply boundary conditions to all molecules in the system.

# Arguments
- `molecules::Vector{<:AbstractOligomer}`: a vector of `AbstractOligomer` objects representing the molecules in the system
- `args::ArgsSmol`: a struct containing the simulation parameters

# Returns
- `nothing`

"""
function apply_boundary!(molecules::Vector{<:AbstractOligomer}, args::ArgsSmol)
    for mol in molecules
        if mol.state == 1 && mol.updated == false
            apply_boundary!(mol, args)
        elseif mol.state == 2
            apply_boundary!(mol, mol.link, args)
        end
    end
    for mol in molecules
        mol.updated = false
    end
end

"""
    apply_boundary!(mol::Monomer, args::ArgsSmol)

Apply boundary conditions to a single monomer in the system.

# Arguments
- `mol::Monomer`: the monomer to apply boundary conditions to
- `args::ArgsSmol`: a struct containing the simulation parameters

# Returns
- `nothing`

"""
function apply_boundary!(mol::Monomer, args::ArgsSmol)
    # Apply boundary conditions to a single molecule
    mol.x, mol.y, mol.z = apply_boundary(mol.x, mol.y, mol.z, args)
    return nothing
end

"""
    apply_boundary!(mol1::Monomer, mol2::Monomer, args::ArgsSmol)

Apply boundary conditions to a dimer using center of mass.

# Note
This may result in molecules being placed outside of the box.

# Arguments
- `mol1::Monomer`: the first monomer in the dimer to apply boundary conditions to
- `mol2::Monomer`: the second monomer in the dimer to apply boundary conditions to
- `args::ArgsSmol`: a struct containing the simulation parameters

# Returns
- `nothing`

"""
function apply_boundary!(mol1::Monomer, mol2::Monomer, args::ArgsSmol)
    # Apply boundary conditions to a dimer using center of mass
    com_x = (mol1.x + mol2.x) / 2
    com_y = (mol1.y + mol2.y) / 2
    com_z = (mol1.z + mol2.z) / 2

    Δx1 = mol1.x - com_x
    Δy1 = mol1.y - com_y
    Δz1 = mol1.z - com_z
    Δx2 = mol2.x - com_x
    Δy2 = mol2.y - com_y
    Δz2 = mol2.z - com_z

    com_x, com_y, com_z = apply_boundary(com_x, com_y, com_z, args)
    mol1.x = com_x + Δx1
    mol1.y = com_y + Δy1
    mol1.z = com_z + Δz1
    mol2.x = com_x + Δx2
    mol2.y = com_y + Δy2
    mol2.z = com_z + Δz2

    return nothing
end

"""
    apply_boundary(x::Float64, y::Float64, z::Float64, args::ArgsSmol)

Apply boundary conditions to a set of coordinates.

# Arguments
- `x::Float64`: the x-coordinate of the point
- `y::Float64`: the y-coordinate of the point
- `z::Float64`: the z-coordinate of the point
- `args::ArgsSmol`: a struct containing the simulation parameters

# Returns
- `x::Float64`: the x-coordinate of the point after applying boundary conditions
- `y::Float64`: the y-coordinate of the point after applying boundary conditions
- `z::Float64`: the z-coordinate of the point after applying boundary conditions

"""
function apply_boundary(x::Float64, y::Float64, z::Float64, args::ArgsSmol)

    box_size_x = args.box_size
    box_size_y = args.box_size
    if args.ndims == 3
        box_size_z = args.box_size
    else
        box_size_z = 0
    end

    if args.boundary == "periodic"
        if x > box_size_x
            x -= box_size_x
        elseif x < 0
            x += box_size_x
        end
        if y > box_size_y
            y -= box_size_y
        elseif y < 0
            y += box_size_y
        end
        if z > box_size_z
            z -= box_size_z
        elseif z < 0
            z += box_size_z
        end
    elseif args.boundary == "reflecting"
        if x > box_size_x
            x = 2 * box_size_x - x
        elseif x < 0
            x = -x
        end
        if y > box_size_y
            y = 2 * box_size_y - y
        elseif y < 0
            y = -y
        end
        if z > box_size_z
            z = 2 * box_size_z - z
        elseif z < 0
            z = -z
        end

    else
        error("Boundary condition not recognized")
    end

    return x, y, z
end


"""
    record_positions!(molecules::Vector{<:AbstractOligomer}, state_history::MoleculeStates, t::Int64)

Record the positions of all molecules in the system at a given time step.

# Arguments
- `molecules::Vector{<:AbstractOligomer}`: a vector of `AbstractOligomer` objects representing the molecules in the system
- `state_history::MoleculeStates`: a `MoleculeStates` object containing the states of the system at each time step
- `t::Int64`: the current time step

# Returns
- `nothing`

"""
function record_positions!(molecules::Vector{<:AbstractOligomer}, state_history::MoleculeStates, t::Int64)
    # Record positions of all molecules
    state_history.States[t] = deepcopy(molecules)
    return nothing
end



"""
    smoluchowski(; kwargs...)

Simulate a system of molecules using the Smoluchowski model.

# Keyword Arguments
- `density::Float64`: the number density of the system (default: 1.0)
- `box_size::Float64`: the size of the simulation box (default: 10.0)
- `diff_monomer::Float64`: the diffusion coefficient of a monomer (default: 0.1)
- `diff_dimer::Float64`: the diffusion coefficient of a dimer (default: 0.05)
- `diff_dimer_rot::Union{Nothing,Float64}`: the rotational diffusion coefficient of a dimer (default: diff_dimer/d_dimer^2)
- `k_off::Float64`: the unbinding rate of a dimer (default: 0.2)
- `r_react::Float64`: the reaction radius for dimerization (default: 0.01)
- `dt::Float64`: the time step used in the simulation (default: 0.01)
- `t_max::Float64`: the maximum simulation time (default: 100.0)
- `ndims::Int64`: the number of dimensions of the system (default: 2)
- `d_dimer::Float64`: the diffusion coefficient of a dimer after unbinding (default: 0.01)
- `boundary::String`: the type of boundary conditions to apply (default: "periodic" or "reflecting")

# Returns
- `state_history::MoleculeStates`: a `MoleculeStates` object containing the states of the system at each time step
- `args::ArgsSmol`: a struct containing the simulation parameters
"""
function smoluchowski(; kwargs...)

    # Parse keyword Arguments
    args = ArgsSmol(; kwargs...)
    # Update orientation
    if isnothing(args.diff_dimer_rot)
        args.diff_dimer_rot = args.diff_dimer/args.d_dimer^2
    end

    # Calculate number of molecules
    n_molecules = round(Int, args.density * args.box_size^args.ndims)
    n_time_steps = round(Int, args.t_max / args.dt)
    # Initialize molecules Simulation States
    state_history = MoleculeStates(args.dt, [Vector{Monomer}(undef, n_molecules) for i in 1:n_time_steps])

    # Initial States
    # These are updated in place during simulation
    molecules = Vector{Monomer}(undef, n_molecules)

    for i in 1:n_molecules
        x = rand(Uniform(0, args.box_size))
        y = rand(Uniform(0, args.box_size))
        if args.ndims == 3
            z = rand(Uniform(0, args.box_size))
        else
            z = 0
        end
        molecules[i] =
            Monomer(x, y, z, 1, nothing, false)
    end

    # Run simulation    
    record_positions!(molecules, state_history, 1)
    for t in 2:n_time_steps
        update_species!(molecules, args)
        update_positions!(molecules, args)
        apply_boundary!(molecules, args)
        record_positions!(molecules, state_history, t)
    end

    return state_history, args
end
