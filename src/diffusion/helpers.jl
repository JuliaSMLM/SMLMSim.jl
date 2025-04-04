"""
    Helpers for diffusion simulation.

This file contains helper functions for the diffusion simulation module.
"""

# Helper functions for geometric calculations
"""
    calc_r(mol1::DiffusingMolecule{<:Emitter2D}, mol2::DiffusingMolecule{<:Emitter2D})

Calculate Euclidean distance between two 2D molecules.

# Arguments
- `mol1::DiffusingMolecule{<:Emitter2D}`: First molecule
- `mol2::DiffusingMolecule{<:Emitter2D}`: Second molecule

# Returns
- `Float64`: Distance between molecules in microns

# Example
```julia
distance = calc_r(molecule1, molecule2)
```
"""
function calc_r(mol1::DiffusingMolecule{<:Emitter2D}, mol2::DiffusingMolecule{<:Emitter2D})
    sqrt((mol1.x - mol2.x)^2 + (mol1.y - mol2.y)^2)
end

"""
    calc_r(mol1::DiffusingMolecule{<:Emitter3D}, mol2::DiffusingMolecule{<:Emitter3D})

Calculate Euclidean distance between two 3D molecules.

# Arguments
- `mol1::DiffusingMolecule{<:Emitter3D}`: First molecule
- `mol2::DiffusingMolecule{<:Emitter3D}`: Second molecule

# Returns
- `Float64`: Distance between molecules in microns
"""
function calc_r(mol1::DiffusingMolecule{<:Emitter3D}, mol2::DiffusingMolecule{<:Emitter3D})
    sqrt((mol1.x - mol2.x)^2 + (mol1.y - mol2.y)^2 + (mol1.z - mol2.z)^2)
end

"""
    calc_ϕ(mol1::DiffusingMolecule, mol2::DiffusingMolecule)

Calculate azimuthal angle between molecules.

# Arguments
- `mol1::DiffusingMolecule`: First molecule
- `mol2::DiffusingMolecule`: Second molecule

# Returns
- `Float64`: Azimuthal angle in radians
"""
function calc_ϕ(mol1::DiffusingMolecule, mol2::DiffusingMolecule)
    atan(mol2.y - mol1.y, mol2.x - mol1.x)
end

"""
    calc_θ(mol1::DiffusingMolecule{<:Emitter2D}, mol2::DiffusingMolecule{<:Emitter2D})

Calculate polar angle between 2D molecules (always π/2 for 2D).

# Arguments
- `mol1::DiffusingMolecule{<:Emitter2D}`: First molecule
- `mol2::DiffusingMolecule{<:Emitter2D}`: Second molecule

# Returns
- `Float64`: Polar angle in radians (always π/2 for 2D)
"""
function calc_θ(mol1::DiffusingMolecule{<:Emitter2D}, mol2::DiffusingMolecule{<:Emitter2D})
    π/2  # Always in xy plane for 2D
end

"""
    calc_θ(mol1::DiffusingMolecule{<:Emitter3D}, mol2::DiffusingMolecule{<:Emitter3D})

Calculate polar angle between 3D molecules.

# Arguments
- `mol1::DiffusingMolecule{<:Emitter3D}`: First molecule
- `mol2::DiffusingMolecule{<:Emitter3D}`: Second molecule

# Returns
- `Float64`: Polar angle in radians
"""
function calc_θ(mol1::DiffusingMolecule{<:Emitter3D}, mol2::DiffusingMolecule{<:Emitter3D})
    r = calc_r(mol1, mol2)
    acos((mol2.z - mol1.z) / r)
end

# Core molecular state changes
"""
    dimerize!(mol1::DiffusingMolecule, mol2::DiffusingMolecule, distance::Real)

Form a dimer between two molecules.

# Arguments
- `mol1::DiffusingMolecule`: First molecule
- `mol2::DiffusingMolecule`: Second molecule
- `distance::Real`: Distance between molecules in the dimer in microns

# Returns
- `Nothing`

# Example
```julia
dimerize!(molecule1, molecule2, 0.05)  # Form dimer with 50nm separation
```

# Note
Molecules are repositioned to maintain the specified distance.
The center of mass position is preserved.
"""
function dimerize!(mol1::DiffusingMolecule, mol2::DiffusingMolecule, distance::Real)
    # Input validation
    if mol1.state == 2 || mol2.state == 2
        return nothing  # Already part of dimers
    end
    
    if distance <= 0
        throw(ArgumentError("Distance must be positive"))
    end
    
    # Update states
    mol1.state = 2
    mol2.state = 2
    mol1.link = mol2.id
    mol2.link = mol1.id
    
    # Calculate center of mass
    com_x = (mol1.x + mol2.x) / 2
    com_y = (mol1.y + mol2.y) / 2
    
    # Calculate orientation
    ϕ = calc_ϕ(mol1, mol2)
    r = distance / 2
    
    # Set new positions based on COM and orientation
    dx = r * cos(ϕ)
    dy = r * sin(ϕ)
    
    # Update positions
    mol1.x = com_x - dx
    mol1.y = com_y - dy
    mol2.x = com_x + dx
    mol2.y = com_y + dy
    
    return nothing
end

"""
    monomerize!(mol::DiffusingMolecule, system::DiffusingMoleculeSystem)

Convert a dimer back to monomers.

# Arguments
- `mol::DiffusingMolecule`: Molecule that is part of a dimer
- `system::DiffusingMoleculeSystem`: System containing all molecules

# Returns
- `Nothing`

# Example
```julia
monomerize!(molecule, system)  # Break dimer into monomers
```

# Note
Both molecules in the dimer are changed to monomer state.
Positions are not changed.
"""
function monomerize!(mol::DiffusingMolecule, system::DiffusingMoleculeSystem)
    if mol.state != 2 || isnothing(mol.link)
        return nothing
    end
    
    # Find linked molecule
    linked_idx = findfirst(m -> m.id == mol.link, system.molecules)
    isnothing(linked_idx) && return nothing
    
    linked_mol = system.molecules[linked_idx]
    
    # Update states
    mol.state = 1
    linked_mol.state = 1
    mol.link = nothing
    linked_mol.link = nothing
    
    return nothing
end