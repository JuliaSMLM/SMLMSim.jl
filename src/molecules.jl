"""
    Molecule

Abstract type for representing photophysical properties of a molecule.

This is the most general type of luminescent or scattering single molecule.
Inherited types will define the properties of specific classes of molecules.
"""
abstract type Molecule end

"""
    GenericFluor <: Molecule

Defines a fluorophore with photophysical properties.

# Fields
- `γ::AbstractFloat`: Photon emission rate in Hz. Default: 1e3
- `q::Array{<:AbstractFloat}`: State transition matrix. Default: q=[1.0]

# Examples
```julia
# Create a fluorophore with default parameters
fluor = GenericFluor()

# Create a fluorophore with custom parameters
fluor = GenericFluor(1e5, [0 10; 1e-1 0])

# Create a fluorophore using keyword arguments
fluor = GenericFluor(; γ=1e5, q=[0 10; 1e-1 0])
```
"""
mutable struct GenericFluor <: Molecule
    γ::AbstractFloat
    q::Array{<:AbstractFloat}
end

function GenericFluor(;
    γ::AbstractFloat=1e5,
    q::Array{<:AbstractFloat}=[1]
)
    return GenericFluor(γ, q)
end