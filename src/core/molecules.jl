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
- `γ::AbstractFloat`: Photon emission rate in Hz. Default: 1e5
- `q::Array{<:AbstractFloat}`: Rate matrix where q[i,j] for i≠j is the transition rate from state i to j,
  and q[i,i] is the negative exit rate from state i. Default: standard 2-state model with
  on->off rate of 50Hz and off->on rate of 1e-2Hz

# Examples
```julia
# Create a fluorophore with default parameters
fluor = GenericFluor()

# Create a fluorophore with custom parameters
fluor = GenericFluor(1e5, [-50.0 50.0; 1e-2 -1e-2])

# Create a fluorophore using keyword arguments
fluor = GenericFluor(; γ=1e5, q=[-10.0 10.0; 1e-1 -1e-1])
```
"""
mutable struct GenericFluor <: Molecule
    γ::AbstractFloat
    q::Array{<:AbstractFloat}
end

function GenericFluor(;
    γ::AbstractFloat=1e5,
    q::Array{<:AbstractFloat}=[-50.0 50.0; 1e-2 -1e-2]
)
    return GenericFluor(γ, q)
end