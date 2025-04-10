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

"""
    Base.show(io::IO, fluor::GenericFluor)

Custom display method for GenericFluor showing basic properties.
"""
function Base.show(io::IO, fluor::GenericFluor)
    n_states = size(fluor.q, 1)
    print(io, "GenericFluor($(n_states) states, γ=$(fluor.γ) Hz)")
end

"""
    Base.show(io::IO, ::MIME"text/plain", fluor::GenericFluor)

Extended display method for GenericFluor in REPL and other text contexts.
"""
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

"""
    Base.show(io::IO, fluor::GenericFluor)

Custom display method for GenericFluor showing basic properties.
"""
function Base.show(io::IO, fluor::GenericFluor)
    n_states = size(fluor.q, 1)
    print(io, "GenericFluor($(n_states) states, γ=$(fluor.γ) Hz)")
end

"""
    Base.show(io::IO, ::MIME"text/plain", fluor::GenericFluor)

Extended display method for GenericFluor in REPL and other text contexts.
"""
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