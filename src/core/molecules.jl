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
# Create a fluorophore with default parameters (using the 2-state keyword constructor)
fluor = GenericFluor()

# Create a fluorophore with custom parameters using the positional constructor
fluor = GenericFluor(1e5, [-50.0 50.0; 1e-2 -1e-2])

# Create a fluorophore using the 2-state keyword constructor
fluor = GenericFluor(; photons=1e5, k_off=10.0, k_on=1e-1)
```
"""
mutable struct GenericFluor <: Molecule
    γ::AbstractFloat
    q::Array{<:AbstractFloat}
end


"""
    GenericFluor(; photons::AbstractFloat=1e5, k_off::AbstractFloat=50.0, k_on::AbstractFloat=1e-2)

Create a simple two-state (on/off) fluorophore with specified parameters.

# Arguments
- `photons::AbstractFloat`: Photon emission rate in Hz
- `k_off::AbstractFloat`: Off-switching rate (on→off) in Hz
- `k_on::AbstractFloat`: On-switching rate (off→on) in Hz

# Details
Creates a fluorophore with a 2-state model and the specified rates.
State 1 is the on (bright) state, and state 2 is the off (dark) state.
The rate matrix is constructed as: q = [-k_off k_off; k_on -k_on]
"""
function GenericFluor(; 
    photons::AbstractFloat=1e5, 
    k_off::AbstractFloat=50.0, 
    k_on::AbstractFloat=1e-2
)
    # Create rate matrix for 2-state system
    q = [-k_off k_off; k_on -k_on]
    
    return GenericFluor(photons, q)
end

function Base.show(io::IO, fluor::GenericFluor)
    n_states = size(fluor.q, 1)
    print(io, "GenericFluor($(n_states) states, γ=$(fluor.γ) Hz)")
end

function Base.show(io::IO, ::MIME"text/plain", fluor::GenericFluor)
    n_states = size(fluor.q, 1)
    println(io, "GenericFluor with $n_states states:")
    println(io, "  Photon emission rate (γ) = $(fluor.γ) Hz")
    println(io, "  Rate matrix (q):")
    
    # Format the rate matrix with aligned columns
    for i in 1:n_states
        print(io, "    ")
        for j in 1:n_states
            val = fluor.q[i, j]
            # Format the rate value with appropriate precision
            val_str = abs(val) < 0.01 ? @sprintf("%.2e", val) : @sprintf("%.3f", val)
            print(io, lpad(val_str, 10))
        end
        if i < n_states
            println(io)
        end
    end
end